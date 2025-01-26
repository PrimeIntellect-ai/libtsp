#include <cassert>
#include <iomanip>
#include <libtsp.h>
#include <random>
#include <unordered_set>
#include <vector>
#include <iostream>

#ifndef WIN32
#define FORCE_INLINE __attribute__((always_inline)) inline
#else
#define FORCE_INLINE inline __forceinline
#endif

#ifndef WIN32
#define PACKED(...) __VA_ARGS__ __attribute__((packed))
#else
#define PACKED(...) __pragma(pack(push, 1)) __VA_ARGS__ __pragma(pack(pop))
#endif

#define TSP_IS_DEBUG

/// Represents an index into the cost array where a particular node is located
typedef std::ptrdiff_t node_cost_idx;

/// Represents an index into the tour array where a particular node is located
typedef std::ptrdiff_t node_tour_idx;

/// Represents a pair of node tour indices
struct node_tour_idx_pair_t {
    node_tour_idx a;
    node_tour_idx b;
};

namespace {
    // define COST_INFINITY
    constexpr cost_t COST_POSITIVE_INFINITY = std::numeric_limits<cost_t>::infinity();

    [[nodiscard]] std::vector<node_cost_idx>
    GenerateRandomTour(const std::unordered_set<node_cost_idx> &remapped_nodes, uint64_t &seed) {
        std::vector<node_cost_idx> tour{};
        tour.reserve(remapped_nodes.size());
        for (const auto node: remapped_nodes) {
            tour.push_back(node);
        }
        std::mt19937_64 rng(seed);
        std::ranges::shuffle(tour, rng);
        seed = rng();
        return tour;
    }

    [[nodiscard]] std::unordered_set<nodeid_t> GetDistinctNodes(const TspInputGraphDescriptor &graph) {
        std::unordered_set<nodeid_t> nodes{};
        for (size_t i = 0; i < graph.num_edges; ++i) {
            nodes.insert(graph.edges[i].from);
            nodes.insert(graph.edges[i].to);
        }
        return nodes;
    }

    [[nodiscard]] size_t GraphCountNodes(const TspInputGraphDescriptor &graph) {
        std::unordered_set<nodeid_t> nodes{};
        for (size_t i = 0; i < graph.num_edges; ++i) {
            nodes.insert(graph.edges[i].from);
            nodes.insert(graph.edges[i].to);
        }
        return nodes.size();
    }

    class CostMatrix;

    void PrintCostMatrix(const CostMatrix &cost_matrix);

    class CostMatrix {
        cost_t *distances;

    public:
        size_t num_nodes;

        CostMatrix() = default;

        CostMatrix(const CostMatrix &) = delete;

        CostMatrix &operator=(const CostMatrix &) = delete;

        CostMatrix(CostMatrix &&other) noexcept: distances(other.distances),
                                                 num_nodes(other.num_nodes) {
            other.distances = nullptr;
            other.num_nodes = 0;
        }

        [[nodiscard]] static std::pair<CostMatrix, std::unordered_map<nodeid_t, node_cost_idx>> CreateCostMatrix(
            const TspInputGraphDescriptor &graph,
            const std::unordered_set<nodeid_t> &distinct_nodes) {
            const size_t num_nodes = GraphCountNodes(graph);

            // build a map of nodeid_t to nodeidx_t
            std::unordered_map<nodeid_t, node_cost_idx> id_to_idx{};
            for (const auto &node: distinct_nodes) {
                id_to_idx[node] = static_cast<node_cost_idx>(id_to_idx.size());
            }

            CostMatrix mat{};
            mat.num_nodes = num_nodes;
            mat.distances = new cost_t[num_nodes * num_nodes];
            std::fill_n(mat.distances, num_nodes * num_nodes, COST_POSITIVE_INFINITY);

            for (size_t i = 0; i < graph.num_edges; ++i) {
                const auto &[cost, from_id, to_id] = graph.edges[i];
                const node_cost_idx from_idx = id_to_idx.at(from_id);
                const node_cost_idx to_idx = id_to_idx.at(to_id);

                mat.distances[from_idx * num_nodes + to_idx] = cost;

                // default to symmetric if user does not specify otherwise
                if (mat.distances[to_idx * num_nodes + from_idx] == COST_POSITIVE_INFINITY) {
                    mat.distances[to_idx * num_nodes + from_idx] = cost;
                }
            }

            // Fill in the diagonal with zeros
            for (size_t i = 0; i < num_nodes; ++i) {
                mat.distances[i * num_nodes + i] = 0;
            }
            return {std::move(mat), std::move(id_to_idx)};
        }

        [[nodiscard]] FORCE_INLINE cost_t get_cost(const node_cost_idx from, const node_cost_idx to) const {
            return distances[from * num_nodes + to];
        }

        ~CostMatrix() {
            delete[] distances;
        }
    };

    /// Prints the cost matrix to the console. Useful for debugging.
    void PrintCostMatrix(const CostMatrix &cost_matrix) {
        // print header
        const auto printHeader = [&cost_matrix] {
            std::cout << "+";
            for (size_t i = 0; i < cost_matrix.num_nodes; ++i) {
                std::cout << "-----"; // allow 4 chars per number
            }
            std::cout << "+" << std::endl;
        };
        printHeader();
        for (node_cost_idx i = 0; i < cost_matrix.num_nodes; ++i) {
            std::cout << "|";
            for (node_cost_idx j = 0; j < cost_matrix.num_nodes; ++j) {
                std::cout << std::setw(4) << cost_matrix.get_cost(i, j) << " ";
            }
            std::cout << "|" << std::endl;
        }
        printHeader();
    }

    PACKED(struct node_pair_t {
        nodeid_t a;
        nodeid_t b;

        bool operator==(const node_pair_t &other) const {
        return a == other.a && b == other.b;
        }

        bool operator!=(const node_pair_t &other) const {
        return !(*this == other);
        }
        });

    static_assert(sizeof(node_pair_t) == 2 * sizeof(nodeid_t));
}

// hash specialization for node_pair_t
template<>
struct [[maybe_unused]] std::hash<node_pair_t> {
    std::size_t operator()(const node_pair_t &pair) const noexcept {
        return std::hash<nodeid_t>{}(pair.a) ^ std::hash<nodeid_t>{}(pair.b);
    }
};

namespace {
    /// Validate that the graph is a valid asymmetric graph.
    /// A valid asymmetric graph must have at least one edge between each pair of nodes.
    /// There may be a second edge connecting the same nodes in the opposite direction with a different cost.
    /// Costs may not be negative.
    /// It also may not contain self-loops.
    TSPStatus ValidateAsymmetricGraph(const TspInputGraphDescriptor &graph) {
        if (graph.num_edges < 2) {
            return TSP_STATUS_ERROR_INVALID_GRAPH;
        }
        std::unordered_set<node_pair_t> edges{};

        // check self-loops, no negative cost, & build set of edges
        for (size_t i = 0; i < graph.num_edges; ++i) {
            const auto &[cost, from, to] = graph.edges[i];
            if (from == to) {
                return TSP_STATUS_ERROR_INVALID_GRAPH;
            }
            if (cost < 0) {
                return TSP_STATUS_ERROR_INVALID_GRAPH;
            }
            const auto &[_, changed] = edges.emplace(from, to);
            if (!changed) {
                return TSP_STATUS_ERROR_INVALID_GRAPH;
            }
        }

        // check duplicate edges
        {
            const std::unordered_set<nodeid_t> distinct_nodes = GetDistinctNodes(graph);
            const std::unordered_set<nodeid_t> &nodes = distinct_nodes;
            for (const auto &a: nodes) {
                for (const auto &b: nodes) {
                    if (a == b) {
                        continue;
                    }
                    const node_pair_t pair_1{a, b};
                    const node_pair_t pair_2{b, a};
                    if (!edges.contains(pair_1) && !edges.contains(pair_2)) {
                        return TSP_STATUS_ERROR_INVALID_GRAPH;
                    }
                }
            }
        }
        return TSP_STATUS_SUCCESS;
    }

    [[nodiscard]] cost_t ComputeTourCost(const std::vector<node_cost_idx> &tour, const CostMatrix &cost_matrix) {
        cost_t cost = 0;
        for (size_t i = 0; i < tour.size(); ++i) {
            const node_cost_idx from_id = tour[i];
            const node_cost_idx to_id = tour[(i + 1) % tour.size()];
            cost += cost_matrix.get_cost(from_id, to_id);
        }
        return cost;
    }

    /// Describes a two-opt move
    struct two_opt_t {
        node_tour_idx i = -1;
        node_tour_idx j = -1;
    };

    /// Describes a three-opt move
    struct three_opt_t {
        node_tour_idx i = -1;
        node_tour_idx j = -1;
        node_tour_idx k = -1;
        bool reverse_states[3] = {false, false, false};
    };

    /// Inplace applies a 2-Opt move to the given tour.
    FORCE_INLINE void ApplyTwoOptMove(std::vector<node_cost_idx> &tour, const two_opt_t &move) {
        node_tour_idx i = move.i;
        node_tour_idx j = move.j;
#ifdef TSP_IS_DEBUG
        if (i > j) {
            std::swap(i, j);
        }
#endif
        std::reverse(tour.begin() + i + 1, tour.begin() + j + 1);
    }

    /// Computes the delta in cost of applying a 2-Opt move to the given tour.
    template<bool asymmetric>
    [[nodiscard]] FORCE_INLINE cost_t ComputeTwoOptMoveDelta(const std::vector<node_cost_idx> &tour,
                                                             const CostMatrix &cost_matrix,
                                                             node_tour_idx i, node_tour_idx j) {
#ifdef TSP_IS_DEBUG
        if (i > j) {
            std::swap(i, j);
        }
#endif

        const size_t n = tour.size();
        const size_t next = (j + 1) % n;

        const node_cost_idx i_idx = tour[i];
        const node_cost_idx ip1_idx = tour[i + 1];
        const node_cost_idx j_idx = tour[j];
        const node_cost_idx jp1_idx = tour[next];

        // Costs of the old edges that will be removed:
        cost_t old_edges = cost_matrix.get_cost(i_idx, ip1_idx) + cost_matrix.get_cost(j_idx, jp1_idx);

        // Costs of the new edges that will be added:
        cost_t new_edges = cost_matrix.get_cost(i_idx, j_idx) + cost_matrix.get_cost(ip1_idx, jp1_idx);

        if constexpr (asymmetric) {
            // For an ASYMMETRIC TSP, also account for all internal edges
            // in the segment i+1 .. j that are reversed:
            // - For k in [i+1..(j-1)], old edge = tour[k] -> tour[k+1]
            //    new edge = tour[k+1] -> tour[k]
            for (node_tour_idx k = i + 1; k < j; ++k) {
                const node_cost_idx k_id = tour[k];
                const node_cost_idx kp1_id = tour[k + 1];
                old_edges += cost_matrix.get_cost(k_id, kp1_id);
                new_edges += cost_matrix.get_cost(kp1_id, k_id);
            }
        }

        // return delta
        return new_edges - old_edges;
    }
}

/// Describes a two-opt tabu entry
struct two_opt_tabu_t {
    nodeid_t a;
    nodeid_t b;

    /// The iteration in which the entry was added to the tabu list
    uint64_t creation_it;

    bool operator==(const two_opt_tabu_t &other) const {
        return a == other.a && b == other.b;
    }

    bool operator!=(const two_opt_tabu_t &other) const {
        return !(*this == other);
    }
};

/// Describes a three-opt tabu entry
struct three_opt_tabu_t {
    nodeid_t a;
    nodeid_t b;
    nodeid_t c;

    /// The iteration in which the entry was added to the tabu list
    uint64_t creation_it;

    bool operator==(const three_opt_tabu_t &other) const {
        return a == other.a && b == other.b && c == other.c;
    }

    bool operator!=(const three_opt_tabu_t &other) const {
        return !(*this == other);
    }
};

// hash specialization for two_opt_t
template<>
struct [[maybe_unused]] std::hash<two_opt_tabu_t> {
    std::size_t operator()(const two_opt_tabu_t &move) const noexcept {
        return std::hash<nodeid_t>{}(move.a) ^ std::hash<nodeid_t>{}(move.b);
    }
};

// hash specialization for three_opt_t
template<>
struct [[maybe_unused]] std::hash<three_opt_tabu_t> {
    std::size_t operator()(const three_opt_tabu_t &move) const noexcept {
        return std::hash<nodeid_t>{}(move.a) ^ std::hash<nodeid_t>{}(move.b) ^ std::hash<nodeid_t>{}(move.c);
    }
};

namespace {
    [[nodiscard]] FORCE_INLINE bool Is2OptTabu(const std::unordered_set<two_opt_tabu_t> &tabu_list,
                                               const nodeid_t a, const nodeid_t b) {
        return tabu_list.contains({a, b});
    }

    [[nodiscard]] FORCE_INLINE bool Is3OptTabu(const std::unordered_set<three_opt_tabu_t> &tabu_list,
                                               const nodeid_t a, const nodeid_t b, const nodeid_t c) {
        return tabu_list.contains({a, b});
    }

    [[nodiscard]] FORCE_INLINE node_tour_idx_pair_t Get3OptSegmentEndpoint(const std::vector<node_cost_idx> &tour,
                                                                           const node_tour_idx seg_start,
                                                                           const node_tour_idx seg_end,
                                                                           const bool is_reversed) {
        return is_reversed
                   ? node_tour_idx_pair_t{tour[seg_end], tour[seg_start]}
                   : node_tour_idx_pair_t{tour[seg_start], tour[seg_end]};
    }

    [[nodiscard]] cost_t ComputeThreeOptMoveDelta(const std::vector<node_cost_idx> &tour,
                                                  const CostMatrix &cost_matrix,
                                                  const node_tour_idx i, const node_tour_idx j, const node_tour_idx k,
                                                  const bool reverse_states[3]) {
        const size_t n = tour.size();
        cost_t old_sum = 0;
        cost_t new_sum = 0;

        // --- Step A: remove E1, E2, E3 (the boundary edges) ---
        old_sum += cost_matrix.get_cost(tour[i], tour[i + 1]);
        old_sum += cost_matrix.get_cost(tour[j], tour[j + 1]);
        old_sum += cost_matrix.get_cost(tour[k], tour[(k + 1) % n]);

        // Remove any arcs inside segments S1, S2 & S3 that if reversed
        if (reverse_states[0]) {
            for (node_tour_idx x = i + 1; x < j; ++x) {
                old_sum += cost_matrix.get_cost(tour[x], tour[x + 1]);
            }
        }
        if (reverse_states[1]) {
            for (node_tour_idx x = j + 1; x < k; ++x) {
                old_sum += cost_matrix.get_cost(tour[x], tour[x + 1]);
            }
        }
        if (reverse_states[2]) {
            for (node_tour_idx x = k + 1; x < i + n; ++x) {
                old_sum += cost_matrix.get_cost(tour[x % n], tour[(x + 1) % n]);
            }
        }

        // --- Step B: add edges for the new arrangement ---

        // compute the index of the respective segment start and end point as if we had
        // actually in-place reversed them into the array.
        const auto &[s1_start, s1_end] = Get3OptSegmentEndpoint(tour, i + 1, j, reverse_states[0]);
        const auto &[s2_start, s2_end] = Get3OptSegmentEndpoint(tour, j + 1, k, reverse_states[1]);
        const auto &[s3_start, s3_end] = Get3OptSegmentEndpoint(tour, k + 1, i, reverse_states[2]);

        // add bridge edges
        new_sum += cost_matrix.get_cost(tour[i], s1_start);
        new_sum += cost_matrix.get_cost(s1_end, s2_start);
        new_sum += cost_matrix.get_cost(s2_end, s3_start);
        new_sum += cost_matrix.get_cost(s3_end, tour[(k + 1) % n]);

        // Add any arcs inside segments S1, S2 & S3 that if reversed
        if (reverse_states[0]) {
            for (node_tour_idx x = i + 1; x < j; ++x) {
                new_sum += cost_matrix.get_cost(tour[x + 1], tour[x]);
            }
        }
        if (reverse_states[1]) {
            for (node_tour_idx x = j + 1; x < k; ++x) {
                new_sum += cost_matrix.get_cost(tour[x + 1], tour[x]);
            }
        }
        if (reverse_states[2]) {
            for (node_tour_idx x = k + 1; x < i + n; ++x) {
                new_sum += cost_matrix.get_cost(tour[(x + 1) % n], tour[x % n]);
            }
        }
        return new_sum - old_sum;
    }

    void ApplyThreeOptMove(std::vector<node_cost_idx> &tour, const three_opt_t &move) {
    }
} // end anonymous namespace

TSP_EXPORT TSPStatus tspSolveAsymmetric(const TspInputGraphDescriptor *graph,
                                        const TspSolverOptionsDescriptor *solver_options,
                                        TspOutputGraphDescriptor *output_descriptor) {
    if (graph == nullptr || output_descriptor == nullptr || solver_options == nullptr) {
        return TSP_STATUS_ERROR_INVALID_ARG;
    }
    if (graph->edges == nullptr || graph->num_edges == 0) {
        return TSP_STATUS_ERROR_INVALID_ARG;
    }
    if (const auto status = ::ValidateAsymmetricGraph(*graph); status != TSP_STATUS_SUCCESS) {
        return status;
    }

    uint64_t seed = solver_options->num_iterations;

    const std::unordered_set<nodeid_t> distinct_nodes = ::GetDistinctNodes(*graph);
    const auto &[cost_matrix, id_to_idx] = CostMatrix::CreateCostMatrix(*graph, distinct_nodes);
    PrintCostMatrix(cost_matrix);

    // remap nodes from id to cost idx
    std::unordered_set<node_cost_idx> remapped_nodes{};
    std::unordered_map<node_cost_idx, nodeid_t> idx_to_id{};
    for (const auto &[id, idx]: id_to_idx) {
        remapped_nodes.insert(idx);
        idx_to_id[idx] = id;
    }

    std::vector<node_cost_idx> current_tour = ::GenerateRandomTour(remapped_nodes, seed);

    cost_t current_cost = ComputeTourCost(current_tour, cost_matrix);

    std::vector<node_cost_idx> best_tour = current_tour;
    cost_t best_cost = current_cost;

    // tabu lists for 2-Opt & 3-Opt moves
    std::unordered_set<two_opt_tabu_t> tabu_list_opt2{};
    std::unordered_set<three_opt_tabu_t> tabu_list_opt3{};

    for (uint64_t iteration = 0; iteration < solver_options->num_iterations; ++iteration) {
        const size_t tour_length = current_tour.size();

        // explore 2-Opt neighbors
        {
            two_opt_t best_move{};
            cost_t best_delta = COST_POSITIVE_INFINITY;

            // Explore 2-Opt neighbors
            for (node_tour_idx i = 0; i < tour_length - 1; i++) {
                for (node_tour_idx j = i + 2; j < tour_length; j++) {
                    const nodeid_t i_id = current_tour[i];
                    const nodeid_t j_id = current_tour[j];

                    const cost_t delta = ::ComputeTwoOptMoveDelta<true>(current_tour, cost_matrix, i, j);
                    const bool is_improvement = delta < 0;
                    const bool better_than_global_best = current_cost + delta < best_cost;
                    if (!::Is2OptTabu(tabu_list_opt2, i_id, j_id) || is_improvement || better_than_global_best) {
                        if (delta < best_delta) {
                            best_delta = delta;
                            best_move = {i, j};
                        }
                    }
                }
            }

            // Update the current tour
            if (best_move.i != -1) {
#ifdef TSP_IS_DEBUG
                assert(best_move.j != -1);
#endif
                ::ApplyTwoOptMove(current_tour, best_move);
                current_cost += best_delta;
                if (current_cost < best_cost) {
                    best_tour = current_tour;
                    best_cost = current_cost;
                }
                tabu_list_opt2.emplace(current_tour[best_move.i], current_tour[best_move.j], iteration);
            }

            // remove all tabu entries where the tenure has expired
            for (auto it = tabu_list_opt2.begin(); it != tabu_list_opt2.end();) {
                if (iteration - it->creation_it >= solver_options->tabu_tenure) {
                    it = tabu_list_opt2.erase(it);
                } else {
                    ++it;
                }
            }
        }

        // explore 3-Opt neighbors
        {
            three_opt_t best_move{};
            cost_t best_delta = COST_POSITIVE_INFINITY;

            // explore 3-Opt neighbors
            for (node_tour_idx i = 0; i < tour_length - 1; i++) {
                size_t i_zero_offs = i == 0 ? 1 : 0;
                for (node_tour_idx j = i + 2; j < tour_length; j++) {
                    for (node_tour_idx k = j + 2; k < tour_length - i_zero_offs; k++) {
                        constexpr static bool reverse_state_permutations[7][3] = {
                            // all permutations for the segment reverse states
                            // 2 ** 3 = 8 permutations (excluding the identity of all false -> 7 permutations to explore)
                            {false, false, true},
                            {false, true, false},
                            {false, true, true},
                            {true, false, false},
                            {true, false, true},
                            {true, true, false},
                            {true, true, true},
                        };
                        for (const auto &segment_reverse_states: reverse_state_permutations) {
                            const cost_t delta = ComputeThreeOptMoveDelta(
                                current_tour, cost_matrix, i, j, k, segment_reverse_states);
                            const bool is_improvement = delta < 0;
                            const bool better_than_global_best = current_cost + delta < best_cost;
                            const nodeid_t i_id = current_tour[i];
                            const nodeid_t j_id = current_tour[j];
                            const nodeid_t k_id = current_tour[k];
                            if (!::Is3OptTabu(tabu_list_opt3, i_id, j_id, k_id) || is_improvement ||
                                better_than_global_best) {
                                if (delta < best_delta) {
                                    best_delta = delta;
                                    best_move = three_opt_t{
                                        .i = i,
                                        .j = j,
                                        .k = k,
                                        .reverse_states = {
                                            segment_reverse_states[0],
                                            segment_reverse_states[1],
                                            segment_reverse_states[2]
                                        }
                                    };
                                }
                            }
                        }
                    }
                }
            }

            // Update the current tour
            if (best_move.i != -1) {
#ifdef TSP_IS_DEBUG
                assert(best_move.j != -1);
                assert(best_move.k != -1);
#endif
                ::ApplyThreeOptMove(current_tour, best_move);
                current_cost += best_delta;
                if (current_cost < best_cost) {
                    best_tour = current_tour;
                    best_cost = current_cost;
                }
                tabu_list_opt3.emplace(current_tour[best_move.i], current_tour[best_move.j], current_tour[best_move.k],
                                       iteration);
            }

            // remove all tabu entries where the tenure has expired
            for (auto it = tabu_list_opt3.begin(); it != tabu_list_opt3.end();) {
                if (iteration - it->creation_it >= solver_options->tabu_tenure) {
                    it = tabu_list_opt3.erase(it);
                } else {
                    ++it;
                }
            }
        }
    }

    // remap back to node ids
    std::vector<nodeid_t> best_tour_ids{};
    best_tour_ids.reserve(best_tour.size());
    for (const auto &idx: best_tour) {
        best_tour_ids.push_back(idx_to_id.at(idx));
    }

    // recompute cost because the constant back and forth additions and subtractions
    // will have introduced floating point errors
    best_cost = ComputeTourCost(best_tour, cost_matrix);

    // Fill output descriptor
    output_descriptor->num_nodes = best_tour_ids.size();
    output_descriptor->solution_cost = best_cost;
    output_descriptor->tour = new nodeid_t[best_tour_ids.size()];
    std::ranges::copy(best_tour_ids, output_descriptor->tour);

    return TSP_STATUS_SUCCESS;
}

TSP_EXPORT TSPStatus tspSolveSymmetric(const TspInputGraphDescriptor *graph,
                                       const TspSolverOptionsDescriptor *solver_options,
                                       TspOutputGraphDescriptor *output_descriptor) {
    // Basic checks
    if (graph == nullptr || output_descriptor == nullptr || solver_options == nullptr) {
        return TSP_STATUS_ERROR_INVALID_ARG;
    }
    if (graph->edges == nullptr || graph->num_edges == 0) {
        return TSP_STATUS_ERROR_INVALID_GRAPH;
    }

    uint64_t seed = solver_options->num_iterations;

    // dummy implementation
    std::vector<node_cost_idx> tour{};

    // Fill output descriptor
    output_descriptor->num_nodes = tour.size();
    output_descriptor->solution_cost = 0;
    output_descriptor->tour = new nodeid_t[tour.size()];
    std::ranges::copy(tour, output_descriptor->tour);

    return TSP_STATUS_SUCCESS;
}
