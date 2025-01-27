#pragma once

#ifndef __cplusplus
#include <stdint.h>
#include <stddef.h>
#else
#include <cstdint>
#include <cstddef>
#endif


#define TSP_EXPORT extern "C"

typedef enum TSPResult {
    TSP_STATUS_SUCCESS = 0, /**< Operation completed successfully. */
    TSP_STATUS_ERROR_INVALID_GRAPH = 1, /**< The input graph is invalid. */
    TSP_STATUS_ERROR_INVALID_ARG = 2, /**< Invalid argument provided. */
    TSP_STATUS_NO_SOLUTION = 3, /**< No solution exists for the given input. */
    TSP_STATUS_OUT_OF_MEMORY = 4, /**< Out of memory. */
    TSP_STATUS_ERROR_INTERNAL = 5 /**< An internal error occurred. */
} TSPStatus;

typedef float cost_t;
typedef uint64_t nodeid_t;

/**
 * Represents a single edge in a TSP graph.
 */
typedef struct TspInputGraphEdge {
    /**
     * The cost associated with the edge. Must be non-negative.
     */
    cost_t cost = 0.0;

    /**
     * A unique 64-bit unsigned integer representing the identity of the node that the edge originates from.
     */
    nodeid_t from = 0;

    /**
     * A unique 64-bit unsigned integer representing the identity of the node that the edge leads to.
     */
    nodeid_t to = 0;
} TspInputGraphEdge;

/**
 * Input descriptor to the TSP solver.
 */
typedef struct TspInputGraphDescriptor {
    /**
     * The edges in the graph
     */
    TspInputGraphEdge *edges{};

    /**
     * The number of elements in the edges array.
     */
    size_t num_edges{};
} TspInputGraphDescriptor;

/**
 * The type of solution returned by the solver
 */
typedef enum TspSolutionType {
    /**
     * The solution is optimal.
     */
    TSP_SOLUTION_TYPE_OPTIMAL = 0,

    /**
     * The solution is approximate.
     */
    TSP_SOLUTION_TYPE_APPROXIMATE = 1,

    /**
     * When requesting an improved solution, whether the solution was actually improved.
     */
    TSP_SOLUTION_TYPE_IMPROVED = 2,

    /**
     * When requesting an improved solution, whether the solution could not be improved.
     */
    TSP_SOLUTION_TYPE_NO_IMPROVEMENT = 3
} TspSolutionType;

typedef struct TspSolutionDescriptor {
    /**
     * The path that represents the solution to the TSP problem. The path is represented as an array of node IDs.
     */
    nodeid_t *tour{};

    /**
     * The number of elements in the path array.
     * Num nodes will be equal to the number of nodes in the graph for a valid solution.
     * An implicit wrap-around edge is assumed from tour[n-1] to tour[0].
     */
    size_t num_nodes{};

    /**
     * The cost of the solution. Must be non-negative.
     * The returned cost will be the sum of the costs of the edges in the path including the
     * final inferred wrap-around edge from tour[n-1] to tour[0].
     */
    cost_t solution_cost{};

    /**
     * The type of the solution.
     */
    TspSolutionType solution_type{};
} TspSolutionDescriptor;

/**
 * A simple enum for choosing the initial-constructive heuristic.
 */
typedef enum TspInitialHeuristic {

    /**
     * Uses a random initial solution.
     */
    TSP_INIT_RANDOM = 0,

    /**
     * Uses the nearest neighbor heuristic to initialize the solution.
     */
    TSP_INIT_NEAREST_NEIGHBOR = 1,

    /**
     * Uses ant colony optimization to initialize the solution.
     */
    TSP_INIT_ANT_COLONY_OPTIMIZATION = 2,

    /**
     * Chooses the strategy that results in the best initial score.
     */
    TSP_INIT_CHOOSE_BEST_INITIAL_SCORE = 3,

    /**
     * Will initialize the solution using a random strategy (nearest neighbor, ant colony optimization, or random with equal probability).
     */
    TSP_INIT_RANDOM_STRATEGY = 4
} TspInitialHeuristic;

/**
 * Solver configuration descriptor.
 */
typedef struct TspSolverConfigurationDescriptor {
    /**
     * Whether to attempt to find an exact solution.
     */
    bool attempt_exact = true;

    /**
     * Upper bound for the number of nodes in the graph for which an exact solution is attempted.
     */
    size_t exact_upper_bound = 16;

    /**
     * The seed for the random number generator
     */
    uint64_t seed{};

    /**
     * The number of iterations to run the local search.
     */
    uint64_t num_iterations = 1000;

    /**
     * The number of iterations that tabu records are retained.
     */
    uint32_t tabu_tenure = 5;

    /**
     * The number of times to restart from a new initial solution.
     * We keep the best solution across all restarts.
     * For each restart, num_iterations are performed.
     * Only meaningful with non-deterministic initialization heuristic.
     */
    uint32_t num_restarts = 1;

    /**
     * Which initial heuristic to use (random, nearest neighbor, or ant colony optimization).
     */
    TspInitialHeuristic initial_heuristic = TSP_INIT_NEAREST_NEIGHBOR;

    /**
     * Number of ants to use in the ant colony optimization heuristic.
     * Only applicable if initial_heuristic is TSP_INIT_ANT_COLONY.
     */
    uint32_t ant_colony_num_samples = 2048;

    /**
     * Whether to perform 3-Opt moves during the local search.
     */
    bool enable_3opt = true;

    /**
     * Whether to perform 4-Opt moves during the local search.
     */
    bool enable_4opt = false;
} TspSolverOptionsDescriptor;

/**
 * Solves the asymmetrical TSP problem for the given graph in an exact manner or heuristically depending on configuration and problem size.
 * @param graph the input graph
 * @param solver_options options to configure the solver
 * @param output_descriptor the output descriptor to write the solution to
 * @return the result status of the operation
 */
TSP_EXPORT TSPStatus tspAsymmetricSolve(const TspInputGraphDescriptor *graph,
                                        const TspSolverOptionsDescriptor *solver_options,
                                        TspSolutionDescriptor *output_descriptor);

/**
 * Improves upon an existing solution provided for an asymmetrical TSP problem in an exact manner or heuristically depending on configuration and problem size.
 * @param graph the input graph
 * @param initial_solution the initial solution to improve upon
 * @param solver_options options to configure the solver
 * @param output_descriptor the output descriptor to write the solution to
 * @return the result status of the operation
 */
TSP_EXPORT TSPStatus tspAsymmetricImproveSolution(const TspInputGraphDescriptor *graph,
                                                  const TspSolutionDescriptor *initial_solution,
                                                  const TspSolverOptionsDescriptor *solver_options,
                                                  TspSolutionDescriptor *output_descriptor);
