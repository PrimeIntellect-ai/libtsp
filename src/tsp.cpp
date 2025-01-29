#include <libtsp.h>
#include <unordered_set>

TSPStatus tspAsymmetricSolveExact(const TspInputGraphDescriptor *graph,
                                  TspSolutionDescriptor *output_descriptor);

TSPStatus tspAsymmetricImproveSolutionHeuristic(const TspInputGraphDescriptor *graph,
                                                const TspSolutionDescriptor *initial_solution,
                                                const TspSolverOptionsDescriptor *solver_options,
                                                TspSolutionDescriptor *output_descriptor);

/// Returns the set of distinct nodes in the graph.
[[nodiscard]] std::unordered_set<nodeid_t> GetDistinctNodes(const TspInputGraphDescriptor &graph) {
    std::unordered_set<nodeid_t> nodes{};
    for (size_t i = 0; i < graph.num_edges; ++i) {
        nodes.insert(graph.edges[i].from);
        nodes.insert(graph.edges[i].to);
    }
    return nodes;
}

struct node_pair_t {
    nodeid_t a;
    nodeid_t b;

    bool operator==(const node_pair_t &other) const {
        return a == other.a && b == other.b;
    }

    bool operator!=(const node_pair_t &other) const {
        return !(*this == other);
    }
};

// hash specialization for node_pair_t
template<>
struct [[maybe_unused]] std::hash<node_pair_t> {
    std::size_t operator()(const node_pair_t &pair) const noexcept {
        return std::hash<nodeid_t>{}(pair.a) ^ std::hash<nodeid_t>{}(pair.b);
    }
};

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
        const auto &[_, changed] = edges.emplace(node_pair_t{from, to});
        if (!changed) {
            return TSP_STATUS_ERROR_INVALID_GRAPH;
        }
    }

    // check that for every pair (a, b), at least one direction is present
    {
        const std::unordered_set<nodeid_t> distinct_nodes = GetDistinctNodes(graph);
        for (const auto &a: distinct_nodes) {
            for (const auto &b: distinct_nodes) {
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

TSPStatus tspAsymmetricSolve(const TspInputGraphDescriptor *graph,
                             const TspSolverOptionsDescriptor *solver_options,
                             TspSolutionDescriptor *output_descriptor) {
    if (graph == nullptr || output_descriptor == nullptr || solver_options == nullptr) {
        return TSP_STATUS_ERROR_INVALID_ARG;
    }
    if (graph->edges == nullptr || graph->num_edges == 0) {
        return TSP_STATUS_ERROR_INVALID_ARG;
    }
    if (const auto status = ::ValidateAsymmetricGraph(*graph); status != TSP_STATUS_SUCCESS) {
        return status;
    }

    const std::unordered_set<nodeid_t> distinct_nodes = GetDistinctNodes(*graph);
    if (solver_options->attempt_exact && distinct_nodes.size() < solver_options->exact_upper_bound) {
        return tspAsymmetricSolveExact(graph, output_descriptor);
    }

    // fall back to heuristic solution
    return tspAsymmetricImproveSolutionHeuristic(graph, nullptr, solver_options, output_descriptor);
}

TSPStatus tspAsymmetricImproveSolution(const TspInputGraphDescriptor *graph,
                                       const TspSolutionDescriptor *initial_solution,
                                       const TspSolverOptionsDescriptor *solver_options,
                                       TspSolutionDescriptor *output_descriptor) {
    if (graph == nullptr || output_descriptor == nullptr || solver_options == nullptr) {
        return TSP_STATUS_ERROR_INVALID_ARG;
    }
    if (graph->edges == nullptr || graph->num_edges == 0) {
        return TSP_STATUS_ERROR_INVALID_ARG;
    }
    if (const auto status = ::ValidateAsymmetricGraph(*graph); status != TSP_STATUS_SUCCESS) {
        return status;
    }

    const std::unordered_set<nodeid_t> distinct_nodes = GetDistinctNodes(*graph);
    if (solver_options->attempt_exact && distinct_nodes.size() < solver_options->exact_upper_bound) {
        return tspAsymmetricSolveExact(graph, output_descriptor);
    }
    return tspAsymmetricImproveSolutionHeuristic(graph, initial_solution, solver_options, output_descriptor);
}

void tspDisposeSolution(const TspSolutionDescriptor *solution) {
    delete[] solution->tour;
}
