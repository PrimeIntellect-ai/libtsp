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
    TSP_STATUS_NO_SOLUTION = 3, /**< No solution found for the given input. */
    TSP_STATUS_OUT_OF_MEMORY = 4, /**< Out of memory. */
    TSP_STATUS_ERROR_INTERNAL = 5 /**< An internal error occurred. */

} TSPStatus;

typedef float cost_t;

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
    uint64_t from = 0;

    /**
     * A unique 64-bit unsigned integer representing the identity of the node that the edge leads to.
     */
    uint64_t to = 0;
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
    uint64_t *path{};

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
 * Solves the TSP problem for the given graph. The solution may not be optimal for large n.
 * @param graph the input graph
 * @param output_descriptor the output descriptor to write the solution to
 * @return the result status of the operation
 */
TSP_EXPORT TSPStatus tsp_solve(const TspInputGraphDescriptor *graph, TspOutputGraphDescriptor *output_descriptor);
