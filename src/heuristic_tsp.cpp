#include "libtsp.h"

#include <memory>
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iostream>

#ifndef WIN32
#define FORCE_INLINE __attribute__((always_inline)) inline
#else
#define FORCE_INLINE inline __forceinline
#endif

/// Represents an index into the cost array where a particular node is located
typedef std::ptrdiff_t node_cost_idx;

/// Represents an index into the tour array where a particular node is located
typedef std::ptrdiff_t node_tour_idx;

/// Represents a pair of node tour indices
struct node_tour_idx_pair_t {
    node_tour_idx a;
    node_tour_idx b;
};

[[nodiscard]] std::unordered_set<nodeid_t> GetDistinctNodes(const TspInputGraphDescriptor &graph);

namespace {
    // define COST_INFINITY
    constexpr cost_t COST_POSITIVE_INFINITY = std::numeric_limits<cost_t>::infinity();

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
        std::unique_ptr<cost_t[]> distances;

    public:
        size_t num_nodes;

        explicit CostMatrix(const size_t num_nodes)
            : distances(std::make_unique<cost_t[]>(num_nodes * num_nodes)),
              num_nodes(num_nodes) {
            std::fill_n(distances.get(), num_nodes * num_nodes, COST_POSITIVE_INFINITY);
        }

        CostMatrix(const CostMatrix &) = delete;

        CostMatrix &operator=(const CostMatrix &) = delete;

        CostMatrix(CostMatrix &&other) noexcept
            : distances(std::move(other.distances)),
              num_nodes(other.num_nodes) {
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

            CostMatrix mat{num_nodes};

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
    };

    /// Prints the cost matrix to the console. Useful for debugging.
    void PrintCostMatrix(const CostMatrix &cost_matrix) {
        // print header
        const auto printHeader = [&cost_matrix] {
            std::cout << "+";
            for (size_t i = 0; i < cost_matrix.num_nodes; ++i) {
                std::cout << "------"; // allow 5 chars per number
            }
            std::cout << "+" << std::endl;
        };
        printHeader();
        for (node_cost_idx i = 0; i < cost_matrix.num_nodes; ++i) {
            std::cout << "|";
            for (node_cost_idx j = 0; j < cost_matrix.num_nodes; ++j) {
                // fix total width to 5 chars
                std::cout << std::fixed << std::setprecision(1) << std::setw(5)
                        << cost_matrix.get_cost(i, j) << " ";
            }
            std::cout << "|" << std::endl;
        }
        printHeader();
    }

    [[nodiscard]] std::vector<node_cost_idx>
    GenerateRandomTour(const std::unordered_set<node_cost_idx> &nodes, uint64_t &seed) {
        std::vector<node_cost_idx> tour{};
        tour.reserve(nodes.size());
        for (const auto node: nodes) {
            tour.push_back(node);
        }
        std::mt19937_64 rng(seed);
        std::ranges::shuffle(tour, rng);
        seed = rng();
        return tour;
    }

    [[nodiscard]] std::vector<node_cost_idx>
    GetNearestNeighborTour(const std::unordered_set<node_cost_idx> &nodes_set,
                           const CostMatrix &cost_matrix, uint64_t &seed) {
        std::vector<node_cost_idx> all_nodes{};
        all_nodes.reserve(nodes_set.size());
        for (const auto &node: nodes_set) {
            all_nodes.push_back(node);
        }

        std::vector<node_cost_idx> tour{};
        tour.reserve(nodes_set.size());

        node_cost_idx starting_point;

        // pick random starting point
        {
            std::mt19937_64 rng(seed);
            std::uniform_int_distribution<size_t> dist(0, all_nodes.size() - 1);
            starting_point = all_nodes[dist(rng)];
            seed = rng();
        }

        tour.push_back(starting_point);

        // build the tour
        while (tour.size() < all_nodes.size()) {
            const node_cost_idx last_node = tour.back();
            cost_t min_cost = COST_POSITIVE_INFINITY;
            node_cost_idx nearest_node = -1;
            for (const auto &node: all_nodes) {
                if (std::ranges::find(tour, node) != tour.end()) {
                    continue;
                }
                const cost_t cost = cost_matrix.get_cost(last_node, node);
                if (cost < min_cost) {
                    min_cost = cost;
                    nearest_node = node;
                }
            }
            assert(nearest_node != -1);
            tour.push_back(nearest_node);
        }
        return tour;
    }

    [[nodiscard]] cost_t ComputeTourCost(const std::vector<node_cost_idx> &tour,
                                         const CostMatrix &cost_matrix) {
        cost_t cost = 0;
        for (size_t i = 0; i < tour.size(); ++i) {
            const node_cost_idx from_id = tour[i];
            const node_cost_idx to_id = tour[(i + 1) % tour.size()];
            cost += cost_matrix.get_cost(from_id, to_id);
        }
        return cost;
    }

    class RewardMatrix {
        std::unique_ptr<float[]> rewards;
        size_t num_nodes;

    public:
        explicit RewardMatrix(const size_t num_nodes)
            : rewards(std::make_unique<float[]>(num_nodes * num_nodes)), num_nodes(num_nodes) {
            std::fill_n(rewards.get(), num_nodes * num_nodes, 0);
        }

        [[nodiscard]] FORCE_INLINE float get_reward(const node_cost_idx from, const node_cost_idx to) const {
            return rewards[from * num_nodes + to];
        }

        FORCE_INLINE void deposit_reward(const node_cost_idx from, const node_cost_idx to, const float reward) {
            rewards[from * num_nodes + to] += reward;
        }
    };

    void RunAntColonySample(std::vector<node_cost_idx> &best_tour, cost_t &best_cost,
                            const std::vector<node_cost_idx> &all_nodes,
                            RewardMatrix &rewards, const CostMatrix &cost_matrix,
                            uint64_t &seed) {
        std::vector<node_cost_idx> tour{};
        tour.reserve(all_nodes.size());

        std::mt19937_64 rng(seed);

        std::unordered_set<node_cost_idx> visited{};

        // choose random starting point
        {
            std::uniform_int_distribution<size_t> dist(0, all_nodes.size() - 1);
            const node_cost_idx starting_point = all_nodes[dist(rng)];
            tour.push_back(starting_point);
            visited.insert(starting_point);
        }

        // sample edges proportional to distance and rewards
        cost_t total_cost = 0;
        for (size_t i = 0; i < all_nodes.size(); i++) {
            const node_cost_idx last_node = tour.back();
            float total_prob = 0;
            for (const auto &node: all_nodes) {
                if (visited.contains(node)) {
                    continue;
                }
                const cost_t cost = cost_matrix.get_cost(last_node, node);
                const float reward = rewards.get_reward(last_node, node);
                total_prob += 1.0f / (cost + reward);
            }

            float choice = std::uniform_real_distribution<float>(0, total_prob)(rng);
            for (const auto &node: all_nodes) {
                if (visited.contains(node)) {
                    continue;
                }
                const cost_t cost = cost_matrix.get_cost(last_node, node);
                const float reward = rewards.get_reward(last_node, node);
                choice -= 1.0f / (cost + reward);
                if (choice <= 0) {
                    tour.push_back(node);
                    visited.insert(node);
                    total_cost += cost;
                    break;
                }
            }
        }

        seed = rng();

        if (total_cost < best_cost) {
            best_cost = total_cost;
            best_tour = tour;

            // deposit rewards
            for (size_t i = 0; i < tour.size(); i++) {
                const node_cost_idx from = tour[i];
                const node_cost_idx to = tour[(i + 1) % tour.size()];
                rewards.deposit_reward(from, to, 1.0f / total_cost);
            }
        }
    }

    [[nodiscard]] std::vector<node_cost_idx>
    AntColonyOptimization(const std::unordered_set<node_cost_idx> &nodes_set,
                          const CostMatrix &cost_matrix,
                          const TspSolverOptionsDescriptor &options_descriptor,
                          uint64_t &seed) {
        RewardMatrix rewards{nodes_set.size()};

        std::vector<node_cost_idx> all_nodes{};
        all_nodes.reserve(nodes_set.size());
        for (const auto &node: nodes_set) {
            all_nodes.push_back(node);
        }

        std::vector<node_cost_idx> best_tour;
        cost_t best_cost = COST_POSITIVE_INFINITY;

        const uint32_t num_samples = std::max(1u, options_descriptor.ant_colony_num_samples);

        for (uint32_t i = 0; i < num_samples; i++) {
            RunAntColonySample(best_tour, best_cost, all_nodes, rewards, cost_matrix, seed);
        }

        return best_tour;
    }
} // end anonymous namespace

namespace {
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

    /// Describes a four-opt move
    struct four_opt_t {
        node_tour_idx w = -1;
        node_tour_idx x = -1;
        node_tour_idx y = -1;
        node_tour_idx z = -1;
        bool reverse_states[4] = {false, false, false, false};
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
        cost_t old_edges = cost_matrix.get_cost(i_idx, ip1_idx)
                           + cost_matrix.get_cost(j_idx, jp1_idx);

        // Costs of the new edges that will be added:
        cost_t new_edges = cost_matrix.get_cost(i_idx, j_idx)
                           + cost_matrix.get_cost(ip1_idx, jp1_idx);

        if constexpr (asymmetric) {
            // For an ASYMMETRIC TSP, also account for edges reversed inside (i+1 .. j)
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

struct four_opt_tabu_t {
    nodeid_t a;
    nodeid_t b;
    nodeid_t c;
    nodeid_t d;

    /// The iteration in which the entry was added to the tabu list
    uint64_t creation_it;

    bool operator==(const four_opt_tabu_t &other) const {
        return (a == other.a && b == other.b && c == other.c && d == other.d);
    }

    bool operator!=(const four_opt_tabu_t &other) const {
        return !(*this == other);
    }
};

// hash specialization for two_opt_t
template<>
struct [[maybe_unused]] std::hash<two_opt_tabu_t> {
    std::size_t operator()(const two_opt_tabu_t &move) const noexcept {
        // simple combination of two node IDs
        auto h1 = std::hash<nodeid_t>{}(move.a);
        auto h2 = std::hash<nodeid_t>{}(move.b);
        // typical combination technique
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
    }
};

// hash specialization for three_opt_t
template<>
struct [[maybe_unused]] std::hash<three_opt_tabu_t> {
    std::size_t operator()(const three_opt_tabu_t &move) const noexcept {
        std::size_t h = 0;
        auto combine = [&](nodeid_t val) {
            h ^= std::hash<nodeid_t>{}(val) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        };
        combine(move.a);
        combine(move.b);
        combine(move.c);
        return h;
    }
};

template<>
struct [[maybe_unused]] std::hash<four_opt_tabu_t> {
    std::size_t operator()(const four_opt_tabu_t &move) const noexcept {
        std::size_t h = 0;
        auto combine = [&](nodeid_t val) {
            h ^= std::hash<nodeid_t>{}(val) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        };
        combine(move.a);
        combine(move.b);
        combine(move.c);
        combine(move.d);
        return h;
    }
};

namespace {
    [[nodiscard]] FORCE_INLINE bool Is2OptTabu(const std::unordered_set<two_opt_tabu_t> &tabu_list,
                                               const nodeid_t a, const nodeid_t b) {
        return tabu_list.contains({a, b});
    }

    [[nodiscard]] FORCE_INLINE bool Is3OptTabu(const std::unordered_set<three_opt_tabu_t> &tabu_list,
                                               const nodeid_t a, const nodeid_t b, const nodeid_t c) {
        return tabu_list.contains({a, b, c});
    }

    [[nodiscard]] FORCE_INLINE bool Is4OptTabu(const std::unordered_set<four_opt_tabu_t> &tabu_list,
                                               const nodeid_t a, const nodeid_t b, const nodeid_t c, const nodeid_t d) {
        return tabu_list.contains({a, b, c, d});
    }

    [[nodiscard]] FORCE_INLINE node_tour_idx_pair_t
    Get3OptSegmentEndpoint(const std::vector<node_cost_idx> &tour,
                           const node_tour_idx seg_start,
                           const node_tour_idx seg_end,
                           const bool is_reversed) {
        return is_reversed
                   ? node_tour_idx_pair_t{tour[seg_end], tour[seg_start]}
                   : node_tour_idx_pair_t{tour[seg_start], tour[seg_end]};
    }

    [[nodiscard]] cost_t ComputeThreeOptMoveDelta(
        const std::vector<node_cost_idx> &tour,
        const CostMatrix &cost_matrix,
        const node_tour_idx i,
        const node_tour_idx j,
        const node_tour_idx k,
        const bool reverse_states[3]) {
        const size_t n = tour.size();

        // 1) Measure the OLD cost
        cost_t old_cost = 0.0;
        for (size_t idx = 0; idx < n; ++idx) {
            const node_cost_idx from_id = tour[idx];
            const node_cost_idx to_id = tour[(idx + 1) % n];
            old_cost += cost_matrix.get_cost(from_id, to_id);
        }

        // 2) Make a scratch copy and APPLY the flips
        static std::vector<node_cost_idx> new_tour{};
        new_tour = tour;

        if (reverse_states[0]) {
            std::reverse(new_tour.begin() + i + 1, new_tour.begin() + j + 1);
        }
        if (reverse_states[1]) {
            std::reverse(new_tour.begin() + j + 1, new_tour.begin() + k + 1);
        }
        if (reverse_states[2]) {
            // wrap-around segment: [k+1 .. end, 0 .. i]
            std::reverse(new_tour.begin() + k + 1, new_tour.end());
            std::reverse(new_tour.begin(), new_tour.begin() + i + 1);
        }

        // 3) Measure the NEW cost
        cost_t new_cost = 0.0;
        for (size_t idx = 0; idx < n; ++idx) {
            const node_cost_idx from_id = new_tour[idx];
            const node_cost_idx to_id = new_tour[(idx + 1) % n];
            new_cost += cost_matrix.get_cost(from_id, to_id);
        }

        // 4) Return delta
        return new_cost - old_cost;
    }

    void ApplyThreeOptMove(std::vector<node_cost_idx> &tour, const three_opt_t &move) {
        const node_tour_idx i = move.i;
        const node_tour_idx j = move.j;
        const node_tour_idx k = move.k;

        if (move.reverse_states[0]) {
            std::reverse(tour.begin() + i + 1, tour.begin() + j + 1);
        }
        if (move.reverse_states[1]) {
            std::reverse(tour.begin() + j + 1, tour.begin() + k + 1);
        }
        if (move.reverse_states[2]) {
            std::reverse(tour.begin() + k + 1, tour.end());
            std::reverse(tour.begin(), tour.begin() + i + 1);
        }
    }

    [[nodiscard]] cost_t ComputeFourOptMoveDelta(
        const std::vector<node_cost_idx> &tour,
        const CostMatrix &cost_matrix,
        const node_tour_idx w,
        const node_tour_idx x,
        const node_tour_idx y,
        const node_tour_idx z,
        const bool reverse_states[4]) {
        const size_t n = tour.size();

        // 1) Old cost
        cost_t old_cost = 0.0;
        for (size_t idx = 0; idx < n; ++idx) {
            const node_cost_idx from_id = tour[idx];
            const node_cost_idx to_id = tour[(idx + 1) % n];
            old_cost += cost_matrix.get_cost(from_id, to_id);
        }

        // 2) Apply flips on a copy
        static std::vector<node_cost_idx> new_tour{};
        new_tour = tour;

        // Segment #0: (w+1 .. x)
        if (reverse_states[0]) {
            std::reverse(new_tour.begin() + w + 1, new_tour.begin() + x + 1);
        }
        // Segment #1: (x+1 .. y)
        if (reverse_states[1]) {
            std::reverse(new_tour.begin() + x + 1, new_tour.begin() + y + 1);
        }
        // Segment #2: (y+1 .. z)
        if (reverse_states[2]) {
            std::reverse(new_tour.begin() + y + 1, new_tour.begin() + z + 1);
        }
        // Segment #3: wrap-around (z+1 .. end) + (0 .. w)
        if (reverse_states[3]) {
            std::reverse(new_tour.begin() + z + 1, new_tour.end());
            std::reverse(new_tour.begin(), new_tour.begin() + w + 1);
        }

        // 3) New cost
        cost_t new_cost = 0.0;
        for (size_t idx = 0; idx < n; ++idx) {
            const node_cost_idx from_id = new_tour[idx];
            const node_cost_idx to_id = new_tour[(idx + 1) % n];
            new_cost += cost_matrix.get_cost(from_id, to_id);
        }

        // 4) Delta
        return new_cost - old_cost;
    }

    void ApplyFourOptMove(std::vector<node_cost_idx> &tour, const four_opt_t &move) {
        const node_tour_idx w = move.w;
        const node_tour_idx x = move.x;
        const node_tour_idx y = move.y;
        const node_tour_idx z = move.z;

        // Segment #0
        if (move.reverse_states[0]) {
            std::reverse(tour.begin() + w + 1, tour.begin() + x + 1);
        }
        // Segment #1
        if (move.reverse_states[1]) {
            std::reverse(tour.begin() + x + 1, tour.begin() + y + 1);
        }
        // Segment #2
        if (move.reverse_states[2]) {
            std::reverse(tour.begin() + y + 1, tour.begin() + z + 1);
        }
        // Segment #3 (wrap-around)
        if (move.reverse_states[3]) {
            std::reverse(tour.begin() + z + 1, tour.end());
            std::reverse(tour.begin(), tour.begin() + w + 1);
        }
    }
} // end anonymous namespace

TSPStatus tspAsymmetricImproveSolutionHeuristic(const TspInputGraphDescriptor *graph,
                                                const TspSolutionDescriptor *initial_solution,
                                                const TspSolverOptionsDescriptor *solver_options,
                                                TspSolutionDescriptor *output_descriptor) {
    uint64_t seed = solver_options->seed;

    const std::unordered_set<nodeid_t> distinct_nodes = ::GetDistinctNodes(*graph);
    const auto &[cost_matrix, id_to_idx] = CostMatrix::CreateCostMatrix(*graph, distinct_nodes);
    // PrintCostMatrix(cost_matrix);

    // remap nodes from id to cost idx
    std::unordered_set<node_cost_idx> remapped_nodes{};
    std::unordered_map<node_cost_idx, nodeid_t> idx_to_id{};
    for (const auto &[id, idx]: id_to_idx) {
        remapped_nodes.insert(idx);
        idx_to_id[idx] = id;
    }

    std::vector<node_cost_idx> best_tour{};
    cost_t best_cost = COST_POSITIVE_INFINITY;
    cost_t initial_cost = COST_POSITIVE_INFINITY;
    if (initial_solution != nullptr) {
        best_cost = initial_cost = initial_solution->solution_cost;
        best_tour.reserve(initial_solution->num_nodes);
        for (size_t i = 0; i < initial_solution->num_nodes; ++i) {
            best_tour.push_back(id_to_idx.at(initial_solution->tour[i]));
        }
    }

    uint32_t num_restarts = std::max(1u, solver_options->num_restarts);
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t restart = 0; restart < num_restarts; restart++) {
        std::vector<node_cost_idx> current_tour;
        cost_t current_cost;
        TspInitialHeuristic chosen_heuristic = solver_options->initial_heuristic;
        if (chosen_heuristic == TSP_INIT_RANDOM_STRATEGY) {
            std::mt19937_64 rng(seed);
            TspInitialHeuristic available_heuristics[] = {
                TSP_INIT_RANDOM,
                TSP_INIT_NEAREST_NEIGHBOR,
                TSP_INIT_ANT_COLONY_OPTIMIZATION,
            };
            std::uniform_int_distribution<size_t> dist(0, 2);
            chosen_heuristic = available_heuristics[dist(rng)];

            seed = rng();
        }
        switch (chosen_heuristic) {
            case TSP_INIT_RANDOM: {
                current_tour = ::GenerateRandomTour(remapped_nodes, seed);
                current_cost = ComputeTourCost(current_tour, cost_matrix);
                break;
            }
            case TSP_INIT_NEAREST_NEIGHBOR: {
                current_tour = ::GetNearestNeighborTour(remapped_nodes, cost_matrix, seed);
                current_cost = ComputeTourCost(current_tour, cost_matrix);
                break;
            }
            case TSP_INIT_ANT_COLONY_OPTIMIZATION: {
                current_tour = ::AntColonyOptimization(remapped_nodes, cost_matrix, *solver_options, seed);
                current_cost = ComputeTourCost(current_tour, cost_matrix);
                break;
            }
            case TSP_INIT_CHOOSE_BEST_INITIAL_SCORE: {
                const auto random_tour = ::GenerateRandomTour(remapped_nodes, seed);
                const auto nearest_neighbor_tour = ::GetNearestNeighborTour(remapped_nodes, cost_matrix, seed);
                const auto ant_colony_tour =
                        ::AntColonyOptimization(remapped_nodes, cost_matrix, *solver_options, seed);

                const cost_t random_cost = ComputeTourCost(random_tour, cost_matrix);
                const cost_t nearest_neighbor_cost = ComputeTourCost(nearest_neighbor_tour, cost_matrix);
                const cost_t ant_colony_cost = ComputeTourCost(ant_colony_tour, cost_matrix);

                if (random_cost < nearest_neighbor_cost && random_cost < ant_colony_cost) {
                    current_tour = random_tour;
                    current_cost = random_cost;
                } else if (nearest_neighbor_cost < random_cost && nearest_neighbor_cost < ant_colony_cost) {
                    current_tour = nearest_neighbor_tour;
                    current_cost = nearest_neighbor_cost;
                } else {
                    current_tour = ant_colony_tour;
                    current_cost = ant_colony_cost;
                }
                break;
            }
            default: {
                return TSP_STATUS_ERROR_INVALID_ARG;
            }
        }

        // The number of times we updated our current best cost via deltas since
        // last force re-computation. We do occasional full cost re-computations
        // to keep floating point drift in check.
        size_t num_dirty_cost_computations = 0;

        if (current_cost < best_cost) {
            best_tour = current_tour;
            best_cost = current_cost;
        }

        // tabu lists for 2-Opt, 3-Opt, and 4-Opt
        std::unordered_set<two_opt_tabu_t> tabu_list_opt2{};
        std::unordered_set<three_opt_tabu_t> tabu_list_opt3{};
        std::unordered_set<four_opt_tabu_t> tabu_list_opt4{};

        for (uint64_t iteration = 0; iteration < solver_options->num_iterations; ++iteration) {
            const size_t tour_length = current_tour.size();

            // explore 2-Opt neighbors
            {
                two_opt_t best_move{};
                cost_t best_delta = COST_POSITIVE_INFINITY;

                for (node_tour_idx i = 0; i < tour_length - 1; i++) {
                    for (node_tour_idx j = i + 2; j < tour_length; j++) {
                        const nodeid_t i_id = current_tour[i];
                        const nodeid_t j_id = current_tour[j];

                        const cost_t delta =
                                ::ComputeTwoOptMoveDelta<true>(current_tour, cost_matrix, i, j);
                        const bool is_improvement = (delta < 0);
                        const bool better_than_global_best = (current_cost + delta < best_cost);

                        if (!::Is2OptTabu(tabu_list_opt2, i_id, j_id)
                            || is_improvement
                            || better_than_global_best) {
                            if (delta < best_delta) {
                                best_delta = delta;
                                best_move = {i, j};
                            }
                        }
                    }
                }

                // Update the current tour if we found a beneficial or permissible 2-Opt move
                if (best_move.i != -1) {
#ifdef TSP_IS_DEBUG
                    assert(best_move.j != -1);
#endif
                    ::ApplyTwoOptMove(current_tour, best_move);
                    current_cost += best_delta;
                    ++num_dirty_cost_computations;

                    // occasionally do a full cost recompute
                    if (num_dirty_cost_computations > 64) {
                        num_dirty_cost_computations = 0;
                        current_cost = ComputeTourCost(current_tour, cost_matrix);
                    }

#ifdef TSP_IS_DEBUG
                    cost_t expected_cost = ComputeTourCost(current_tour, cost_matrix);
                    assert(fabs(current_cost - expected_cost) < 1e-3);
#endif

                    if (current_cost < best_cost) {
                        best_tour = current_tour;
                        best_cost = current_cost;
                    }
                    tabu_list_opt2.emplace(
                        current_tour[best_move.i],
                        current_tour[best_move.j],
                        iteration
                    );
                }

                // remove expired 2-Opt tabu
                for (auto it = tabu_list_opt2.begin(); it != tabu_list_opt2.end();) {
                    if (iteration - it->creation_it >= solver_options->tabu_tenure) {
                        it = tabu_list_opt2.erase(it);
                    } else {
                        ++it;
                    }
                }
            }

            // explore 3-Opt neighbors
            if (solver_options->enable_3opt) {
                three_opt_t best_move{};
                cost_t best_delta = COST_POSITIVE_INFINITY;

                for (node_tour_idx i = 0; i < tour_length - 1; i++) {
                    size_t i_zero_offs = (i == 0) ? 1 : 0;
                    for (node_tour_idx j = i + 2; j < tour_length; j++) {
                        for (node_tour_idx k = j + 2; k < tour_length - i_zero_offs; k++) {
                            // 2^3 = 8 permutations minus identity => 7
                            constexpr static bool reverse_state_permutations[7][3] = {
                                {false, false, true}, {false, true, false}, {false, true, true},
                                {true, false, false}, {true, false, true}, {true, true, false},
                                {true, true, true},
                            };
                            for (const auto &segment_reverse_states: reverse_state_permutations) {
                                const cost_t delta = ComputeThreeOptMoveDelta(
                                    current_tour, cost_matrix, i, j, k, segment_reverse_states);

                                const bool is_improvement = (delta < 0);
                                const bool better_than_global_best =
                                        (current_cost + delta < best_cost);

                                const nodeid_t i_id = current_tour[i];
                                const nodeid_t j_id = current_tour[j];
                                const nodeid_t k_id = current_tour[k];

                                if (!::Is3OptTabu(tabu_list_opt3, i_id, j_id, k_id)
                                    || is_improvement
                                    || better_than_global_best) {
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

                // Update the current tour if 3-Opt move found
                if (best_move.i != -1) {
#ifdef TSP_IS_DEBUG
                    assert(best_move.j != -1);
                    assert(best_move.k != -1);
#endif
                    ::ApplyThreeOptMove(current_tour, best_move);
                    current_cost += best_delta;
                    ++num_dirty_cost_computations;

                    if (num_dirty_cost_computations > 64) {
                        num_dirty_cost_computations = 0;
                        current_cost = ComputeTourCost(current_tour, cost_matrix);
                    }

#ifdef TSP_IS_DEBUG
                    cost_t expected_cost = ComputeTourCost(current_tour, cost_matrix);
                    assert(fabs(current_cost - expected_cost) < 1e-3);
#endif

                    if (current_cost < best_cost) {
                        best_tour = current_tour;
                        best_cost = current_cost;
                    }

                    tabu_list_opt3.emplace(
                        current_tour[best_move.i],
                        current_tour[best_move.j],
                        current_tour[best_move.k],
                        iteration
                    );
                }

                // remove expired 3-Opt tabu
                for (auto it = tabu_list_opt3.begin(); it != tabu_list_opt3.end();) {
                    if (iteration - it->creation_it >= solver_options->tabu_tenure) {
                        it = tabu_list_opt3.erase(it);
                    } else {
                        ++it;
                    }
                }
            }

            // 4-Opt exploration
            if (solver_options->enable_4opt) {
                four_opt_t best_move{};
                cost_t best_delta = COST_POSITIVE_INFINITY;

                for (node_tour_idx w = 0; w < tour_length - 1; w++) {
                    // Similar offset logic if needed:
                    size_t w_zero_offs = (w == 0) ? 1 : 0;
                    for (node_tour_idx x = w + 2; x < tour_length; x++) {
                        for (node_tour_idx y = x + 2; y < tour_length; y++) {
                            for (node_tour_idx z = y + 2; z < tour_length - w_zero_offs; z++) {
                                // We skip the identity (false,false,false,false) => 2^4 - 1 = 15 combos
                                constexpr static bool reverse_state_permutations[15][4] = {
                                    {false, false, false, true},
                                    {false, false, true, false},
                                    {false, false, true, true},
                                    {false, true, false, false},
                                    {false, true, false, true},
                                    {false, true, true, false},
                                    {false, true, true, true},
                                    {true, false, false, false},
                                    {true, false, false, true},
                                    {true, false, true, false},
                                    {true, false, true, true},
                                    {true, true, false, false},
                                    {true, true, false, true},
                                    {true, true, true, false},
                                    {true, true, true, true},
                                };
                                for (const auto &rev_state: reverse_state_permutations) {
                                    const cost_t delta = ComputeFourOptMoveDelta(
                                        current_tour, cost_matrix, w, x, y, z, rev_state);

                                    const bool is_improvement = (delta < 0);
                                    const bool better_than_global_best =
                                            (current_cost + delta < best_cost);

                                    const nodeid_t w_id = current_tour[w];
                                    const nodeid_t x_id = current_tour[x];
                                    const nodeid_t y_id = current_tour[y];
                                    const nodeid_t z_id = current_tour[z];

                                    if (!::Is4OptTabu(tabu_list_opt4, w_id, x_id, y_id, z_id)
                                        || is_improvement
                                        || better_than_global_best) {
                                        if (delta < best_delta) {
                                            best_delta = delta;
                                            best_move = four_opt_t{
                                                .w = w,
                                                .x = x,
                                                .y = y,
                                                .z = z,
                                                .reverse_states = {
                                                    rev_state[0], rev_state[1],
                                                    rev_state[2], rev_state[3]
                                                }
                                            };
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // Update if found a 4-Opt move
                if (best_move.w != -1) {
                    ::ApplyFourOptMove(current_tour, best_move);
                    current_cost += best_delta;
                    ++num_dirty_cost_computations;

                    if (num_dirty_cost_computations > 64) {
                        num_dirty_cost_computations = 0;
                        current_cost = ComputeTourCost(current_tour, cost_matrix);
                    }

#ifdef TSP_IS_DEBUG
                    cost_t expected_cost = ComputeTourCost(current_tour, cost_matrix);
                    assert(fabs(current_cost - expected_cost) < 1e-3);
#endif

                    if (current_cost < best_cost) {
                        best_tour = current_tour;
                        best_cost = current_cost;
                    }

                    // add to 4-Opt tabu
                    tabu_list_opt4.emplace(
                        current_tour[best_move.w],
                        current_tour[best_move.x],
                        current_tour[best_move.y],
                        current_tour[best_move.z],
                        iteration
                    );
                }

                // remove expired 4-Opt tabu
                for (auto it = tabu_list_opt4.begin(); it != tabu_list_opt4.end();) {
                    if (iteration - it->creation_it >= solver_options->tabu_tenure) {
                        it = tabu_list_opt4.erase(it);
                    } else {
                        ++it;
                    }
                }
            }

            if (solver_options->time_limit_ms != UINT64_MAX) {
                auto now = std::chrono::high_resolution_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
                if (elapsed >= solver_options->time_limit_ms) {
                    break;
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

    // final cost recompute for safety
    best_cost = ComputeTourCost(best_tour, cost_matrix);

    // Fill output descriptor
    output_descriptor->num_nodes = best_tour_ids.size();
    output_descriptor->solution_cost = best_cost;
    if (initial_solution == nullptr) {
        output_descriptor->solution_type = TSP_SOLUTION_TYPE_APPROXIMATE;
    } else {
        if (best_cost < initial_cost) {
            output_descriptor->solution_type = TSP_SOLUTION_TYPE_IMPROVED;
        } else {
            output_descriptor->solution_type = TSP_SOLUTION_TYPE_NO_IMPROVEMENT;
        }
    }
    output_descriptor->tour = new nodeid_t[best_tour_ids.size()];
    std::ranges::copy(best_tour_ids, output_descriptor->tour);

    return TSP_STATUS_SUCCESS;
}
