#include "libtsp.h"

#include <vector>
#include <unordered_map>
#include <limits>
#include <cstring>


static constexpr cost_t COST_POSITIVE_INFINITY = std::numeric_limits<cost_t>::infinity();

namespace {
    /**
     * Maps each unique node ID to a continuous index [0..(N-1)].
     * Returns a vector of the unique node IDs in the order they were assigned indices.
     */
    std::vector<nodeid_t> MapNodeIdsToIndices(
        const TspInputGraphDescriptor *graph,
        std::unordered_map<nodeid_t, size_t> &nodeIndexMap
    ) {
        std::vector<nodeid_t> uniqueNodes;
        uniqueNodes.reserve(graph->num_edges * 2); // rough upper bound

        // Collect unique node IDs
        for (size_t i = 0; i < graph->num_edges; i++) {
            nodeid_t fromNode = graph->edges[i].from;
            nodeid_t toNode = graph->edges[i].to;

            // If encountering a new node, give it an index
            if (!nodeIndexMap.contains(fromNode)) {
                const size_t idx = nodeIndexMap.size();
                nodeIndexMap[fromNode] = idx;
                uniqueNodes.push_back(fromNode);
            }
            if (!nodeIndexMap.contains(toNode)) {
                const size_t idx = nodeIndexMap.size();
                nodeIndexMap[toNode] = idx;
                uniqueNodes.push_back(toNode);
            }
        }

        return uniqueNodes;
    }

    /**
     * Constructs an NxN cost matrix, where N is the number of unique nodes.
     * costMatrix[i][j] = cost of edge from node i to node j, or INF if no edge exists.
     */
    std::vector<cost_t> buildCostMatrix(
        const TspInputGraphDescriptor *graph,
        const std::unordered_map<nodeid_t, size_t> &nodeIndexMap,
        const size_t n
    ) {
        // Initialize all costs to INF
        std::vector costMatrix(n * n, COST_POSITIVE_INFINITY);

        // Fill in the actual costs
        // If there are multiple edges between the same pair (from->to),
        // we'll use the minimal cost among them. You can adjust as needed.
        for (size_t i = 0; i < graph->num_edges; i++) {
            const auto &[cost, from, to] = graph->edges[i];
            if (cost < 0.0f) {
                // We do not handle negative costs in this basic implementation;
                // this will be caught outside this function as invalid graph
                continue;
            }
            const size_t fromIndex = nodeIndexMap.at(from);
            const size_t toIndex = nodeIndexMap.at(to);
            costMatrix[fromIndex * n + toIndex] =
                    std::min(costMatrix[fromIndex * n + toIndex], cost);
        }

        return costMatrix;
    }

    /**
     * Performs the Held-Karp dynamic programming algorithm for the ATSP.
     *
     * We choose the first node (index 0) as the "start" node for simplicity.
     * The final route is forced to return back to this start node.
     *
     * DP state: dp[mask][i] = minimal cost to start at node 0, visit all nodes in mask exactly once,
     * ending at node i (where i is in mask).
     *
     * mask is a bitmask representing which nodes have been visited.
     * parent[mask][i] is used to reconstruct the path.
     *
     * N should be the total number of distinct nodes.
     */
    TSPStatus HeldKarpSolve(
        const size_t n,
        const std::vector<cost_t> &costMatrix,
        std::vector<size_t> &bestPath,
        cost_t &bestCost
    ) {
        if (n == 0) {
            // Trivial case: no nodes
            bestCost = 0.0f;
            return TSP_STATUS_NO_SOLUTION;
        }
        if (n == 1) {
            // Trivial case: single node, the tour is that node
            bestPath.clear();
            bestPath.push_back(0);
            bestPath.push_back(0); // return to self for a cycle
            bestCost = 0.0f;
            return TSP_STATUS_SUCCESS;
        }

        // Number of possible subsets is 1 << N
        const size_t subsetCount = (1ULL << n);

        // dp[mask][i]: the minimal cost to start at node 0, visit all nodes in mask,
        // ending at node i. We'll store these in a 2D array with dimension subsetCount x N.
        std::unique_ptr<cost_t[]> dp(new(std::nothrow) cost_t[subsetCount * n]);
        if (!dp) {
            return TSP_STATUS_OUT_OF_MEMORY;
        }
        // parent[mask][i] to reconstruct path: which predecessor node j led to i
        std::unique_ptr<int[]> parent(new(std::nothrow) int[subsetCount * n]);
        if (!parent) {
            return TSP_STATUS_OUT_OF_MEMORY;
        }

        // Initialize dp array to INF
        for (size_t i = 0; i < subsetCount * n; i++) {
            dp[i] = COST_POSITIVE_INFINITY;
            parent[i] = -1;
        }

        // We define a start node to be 0 (the first in our indexing)
        constexpr size_t start = 0;
        constexpr size_t startMask = (1ULL << start);

        // Base case: dp[(1 << start)][start] = 0
        dp[startMask * n + start] = 0.0f;

        // Build up the dp table
        for (size_t mask = 0; mask < subsetCount; mask++) {
            // If the start node is not in this subset, skip
            if ((mask & startMask) == 0) {
                continue;
            }

            // For each possible 'last node' i in this subset
            for (size_t i = 0; i < n; i++) {
                if ((mask & (1ULL << i)) == 0) {
                    continue; // i not in subset
                }

                const cost_t currentCost = dp[mask * n + i];
                if (currentCost == COST_POSITIVE_INFINITY) {
                    continue; // not a valid state
                }

                // Try to extend the path by going to a node j not in 'mask'
                for (size_t j = 0; j < n; j++) {
                    if ((mask & (1ULL << j)) != 0) {
                        continue; // j already in subset
                    }
                    // cost i->j
                    const cost_t edgeCost = costMatrix[i * n + j];
                    if (edgeCost == COST_POSITIVE_INFINITY) {
                        continue; // no edge i->j
                    }
                    const size_t nextMask = mask | (1ULL << j);
                    const cost_t newCost = currentCost + edgeCost;
                    cost_t &dpRef = dp[nextMask * n + j];
                    if (newCost < dpRef) {
                        dpRef = newCost;
                        parent[nextMask * n + j] = static_cast<int>(i);
                    }
                }
            }
        }

        // To get a complete tour, we must return to the start node (0).
        // We consider dp[(1<<N) - 1][i] + cost(i->start).
        const size_t allVisitedMask = (1ULL << n) - 1;
        bestCost = COST_POSITIVE_INFINITY;
        int bestLastNode = -1;

        for (size_t i = 0; i < n; i++) {
            if (i == start) {
                continue; // we skip if i==start, but it's not necessarily invalid to check
            }
            const cost_t routeCost = dp[allVisitedMask * n + i];
            if (routeCost >= COST_POSITIVE_INFINITY) {
                continue; // no route ending in i visiting all nodes
            }
            const cost_t returnEdge = costMatrix[i * n + start];
            if (returnEdge == COST_POSITIVE_INFINITY) {
                continue; // can't return to start from i
            }
            const cost_t totalCost = routeCost + returnEdge;
            if (totalCost < bestCost) {
                bestCost = totalCost;
                bestLastNode = static_cast<int>(i);
            }
        }

        if (bestLastNode < 0 || bestCost >= COST_POSITIVE_INFINITY) {
            // No solution found that visits all nodes and returns to start
            return TSP_STATUS_NO_SOLUTION;
        }

        // Reconstruct path
        bestPath.clear();
        bestPath.resize(n, 0);

        // We'll fill bestPath from the end to the start
        size_t curMask = allVisitedMask;
        int curNode = bestLastNode;

        for (int pos = static_cast<int>(n) - 1; pos >= 1; pos--) {
            bestPath[pos] = static_cast<size_t>(curNode);
            const int prevNode = parent[curMask * n + curNode];
            // Remove curNode from curMask
            curMask = curMask & ~(1ULL << curNode);
            curNode = prevNode;
        }

        // The start node is always 0
        bestPath[0] = start;

        return TSP_STATUS_SUCCESS;
    }
}

TSPStatus tspAsymmetricSolveExact(const TspInputGraphDescriptor *graph, TspSolutionDescriptor *output_descriptor) {
    // Basic argument checking
    if (!graph || !output_descriptor) {
        return TSP_STATUS_ERROR_INVALID_ARG;
    }
    if (!graph->edges || graph->num_edges == 0) {
        // No edges means we cannot form a valid tour
        return TSP_STATUS_ERROR_INVALID_GRAPH;
    }

    // Collect all unique node IDs and map them to [0..N-1]
    std::unordered_map<nodeid_t, size_t> nodeIndexMap;
    const std::vector<nodeid_t> indexToNode = MapNodeIdsToIndices(graph, nodeIndexMap);
    const size_t N = indexToNode.size();

    // Check for negative costs (invalid) and build cost matrix
    for (size_t i = 0; i < graph->num_edges; i++) {
        if (graph->edges[i].cost < 0.0f) {
            return TSP_STATUS_ERROR_INVALID_GRAPH;
        }
    }
    const std::vector<cost_t> costMatrix = buildCostMatrix(graph, nodeIndexMap, N);

    // Now run Held-Karp
    std::vector<size_t> bestPathIndices; // will hold node indices in [0..N-1]
    cost_t bestCost = 0.0f;

    const TSPStatus status = HeldKarpSolve(N, costMatrix, bestPathIndices, bestCost);
    if (status != TSP_STATUS_SUCCESS) {
        return status;
    }

    // Convert bestPathIndices (which are indices) back to node IDs
    // bestPathIndices.size() should be N+1
    std::vector<nodeid_t> bestPathNodeIds(bestPathIndices.size());
    for (size_t i = 0; i < bestPathIndices.size(); i++) {
        const size_t idx = bestPathIndices[i];
        bestPathNodeIds[i] = indexToNode[idx];
    }

    // Prepare output descriptor
    // We allocate memory for the tour. The caller is responsible for freeing it
    // as per typical C-API conventions, or you can adjust to your memory model.
    output_descriptor->num_nodes = bestPathNodeIds.size();
    output_descriptor->tour = new(std::nothrow) nodeid_t[bestPathNodeIds.size()];
    if (!output_descriptor->tour) {
        return TSP_STATUS_OUT_OF_MEMORY;
    }

    // Copy the path
    std::memcpy(output_descriptor->tour,
                bestPathNodeIds.data(),
                bestPathNodeIds.size() * sizeof(nodeid_t));

    // Set cost and solution type
    output_descriptor->solution_cost = bestCost;
    output_descriptor->solution_type = TSP_SOLUTION_TYPE_OPTIMAL;

    return TSP_STATUS_SUCCESS;
}
