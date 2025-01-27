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
    TSP_STATUS_OUT_OF_MEMORY = 3, /**< Out of memory. */
    TSP_STATUS_ERROR_INTERNAL = 4 /**< An internal error occurred. */
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

typedef struct TspOutputGraphDescriptor {
    /**
     * The path that represents the solution to the TSP problem. The path is represented as an array of node IDs.
     */
    nodeid_t *tour{};

    /**
     * The number of elements in the path array.
     */
    size_t num_nodes{};

    /**
     * The cost of the solution. Must be non-negative.
     */
    cost_t solution_cost{};
} TspOutputGraphDescriptor;

/**
 * A simple enum for choosing the initial-constructive heuristic.
 */
typedef enum TspInitialHeuristic {
    TSP_INIT_RANDOM = 0,
    TSP_INIT_NEAREST_NEIGHBOR = 1
} TspInitialHeuristic;

/**
 * Solver configuration descriptor.
 */
typedef struct TspSolverConfigurationDescriptor {
    /**
     * The seed for the random number generator (if needed).
     */
    uint64_t seed{};

    /**
     * The number of iterations to run the local search.
     */
    uint64_t num_iterations = 1000;

    /**
     * The number of iterations that tabu records are retained.
     */
    uint32_t tabu_tenure = 10;

    /**
     * The number of times to restart from a new initial solution.
     * We keep the best solution across all restarts.
     */
    uint32_t num_restarts = 1;

    /**
     * Which initial heuristic to use (random or nearest-neighbor).
     */
    TspInitialHeuristic initial_heuristic = TSP_INIT_NEAREST_NEIGHBOR;

    /**
     * Whether to run 3-Opt after 2-Opt each iteration.
     * If false, only 2-Opt steps are performed.
     */
    bool enable_3opt = true;
} TspSolverOptionsDescriptor;

/**
 * Solves the asymmetrical TSP problem for the given graph.
 * @param graph the input graph
 * @param solver_options options to configure the solver
 * @param output_descriptor the output descriptor to write the solution to
 * @return the result status of the operation
 */
TSP_EXPORT TSPStatus tspSolveAsymmetric(const TspInputGraphDescriptor *graph,
                                        const TspSolverOptionsDescriptor *solver_options,
                                        TspOutputGraphDescriptor *output_descriptor);

/**
 * Solves the symmetrical TSP problem for the given graph.
 * @param graph the input graph
 * @param solver_options options to configure the solver
 * @param output_descriptor the output descriptor to write the solution to
 * @return the result status of the operation
 */
TSP_EXPORT TSPStatus tspSolveSymmetric(const TspInputGraphDescriptor *graph,
                                       const TspSolverOptionsDescriptor *solver_options,
                                       TspOutputGraphDescriptor *output_descriptor);
