#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <memory>
#include <new>
#include <vector>
#include <iostream>

// If using OpenMP, include <omp.h> and compile with -fopenmp or equivalent.
#ifdef _OPENMP
#include <omp.h>
#endif

#include <libtsp.h> // Adjust if needed

static constexpr cost_t INF = std::numeric_limits<cost_t>::infinity();

namespace {

    struct NodeMapper {
        std::vector<uint64_t> index_to_id;
        std::map<uint64_t, int> id_to_index;

        void build(const TspInputGraphDescriptor *graph) {
            index_to_id.reserve(graph->num_edges * 2ULL);
            for (size_t i = 0; i < graph->num_edges; ++i) {
                auto from_id = graph->edges[i].from;
                auto to_id   = graph->edges[i].to;
                if (!id_to_index.contains(from_id)) {
                    const int idx = static_cast<int>(index_to_id.size());
                    index_to_id.push_back(from_id);
                    id_to_index[from_id] = idx;
                }
                if (!id_to_index.contains(to_id)) {
                    const int idx = static_cast<int>(index_to_id.size());
                    index_to_id.push_back(to_id);
                    id_to_index[to_id] = idx;
                }
            }
        }
        [[nodiscard]] size_t size() const {
            return index_to_id.size();
        }
    };

    // Fast popcount for 32-bit.
    int popcount(unsigned x) {
#ifdef __GNUC__
        return __builtin_popcount(x);
#else
        int c = 0;
        while (x) {
            x &= (x - 1);
            c++;
        }
        return c;
#endif
    }

    // Reconstruct path from DP parent info.
    void reconstructPath(
        const int *parent,
        const int n,
        const int fullMask,
        const int endNode,
        std::vector<int> &pathOut)
    {
        pathOut.clear();
        int curMask = fullMask;
        int curNode = endNode;
        while (true) {
            pathOut.push_back(curNode);
            int p = parent[static_cast<size_t>(curMask) * n + curNode];
            if (p < 0) {
                // means we came from node 0
                break;
            }
            curMask ^= (1 << (curNode - 1));
            curNode = p;
        }
        std::ranges::reverse(pathOut);
    }

    struct IntCostPair {
        uint32_t a;
        cost_t b;
    };
    static_assert(sizeof(IntCostPair) == 8);

} // end anonymous namespace


TSP_EXPORT TSPStatus tsp_solve(const TspInputGraphDescriptor *graph,
                               TspOutputGraphDescriptor *output_descriptor)
{
    // Basic checks
    if (!graph || !output_descriptor) {
        return TSP_STATUS_ERROR_INVALID_ARG;
    }
    if (!graph->edges || graph->num_edges == 0) {
        return TSP_STATUS_ERROR_INVALID_ARG;
    }

    NodeMapper mapper;
    mapper.build(graph);
    const int n = static_cast<int>(mapper.size());
    if (n == 0) {
        return TSP_STATUS_NO_SOLUTION;
    }
    if (n == 1) {
        // Only one node => cost=0 path
        output_descriptor->num_nodes = 1;
        output_descriptor->path = new(std::nothrow) uint64_t[1];
        if (!output_descriptor->path) {
            return TSP_STATUS_OUT_OF_MEMORY;
        }
        output_descriptor->path[0] = mapper.index_to_id[0];
        output_descriptor->solution_cost = 0.0;
        return TSP_STATUS_SUCCESS;
    }

    // Build full adjacency matrix
    std::vector<cost_t> adjacency(n * n, INF);
    adjacency.shrink_to_fit();

    for (size_t i = 0; i < graph->num_edges; i++) {
        cost_t c = graph->edges[i].cost;
        if (c < 0.0) {
            return TSP_STATUS_ERROR_INVALID_ARG; // negative cost invalid
        }
        int u = mapper.id_to_index[graph->edges[i].from];
        int v = mapper.id_to_index[graph->edges[i].to];
        if (adjacency[static_cast<size_t>(u) * n + v] != INF) {
            // duplicate edge
            return TSP_STATUS_ERROR_INVALID_GRAPH;
        }
        adjacency[static_cast<size_t>(u) * n + v] = c;
    }

    // Handle n=2 quickly
    if (n == 2) {
        cost_t c1 = adjacency[0 * n + 1];
        cost_t c2 = adjacency[1 * n + 0];
        if (c1 >= INF || c2 >= INF) {
            return TSP_STATUS_NO_SOLUTION;
        }
        output_descriptor->num_nodes = 3;
        auto *sol = new(std::nothrow) uint64_t[3];
        if (!sol) return TSP_STATUS_OUT_OF_MEMORY;
        sol[0] = mapper.index_to_id[0];
        sol[1] = mapper.index_to_id[1];
        sol[2] = mapper.index_to_id[0];
        output_descriptor->path = sol;
        output_descriptor->solution_cost = c1 + c2;
        return TSP_STATUS_SUCCESS;
    }

    // Edge-dominance elimination (O(n^3)), optional but often helpful
    for (int u = 0; u < n; u++) {
        for (int v = 0; v < n; v++) {
            cost_t uv = adjacency[u * n + v];
            if (u != v && uv < INF) {
                for (int w = 0; w < n; w++) {
                    if (w == u || w == v) continue;
                    cost_t uw = adjacency[u * n + w];
                    cost_t wv = adjacency[w * n + v];
                    if (uw < INF && wv < INF) {
                        if (uv > uw + wv) {
                            adjacency[u * n + v] = INF;
                            break;
                        }
                    }
                }
            }
        }
    }

    // Build sparse adjacency lists
    std::vector<std::vector<IntCostPair>> neighbors(n);
    neighbors.shrink_to_fit();
    for (int u = 0; u < n; u++) {
        std::vector<IntCostPair> row;
        row.reserve(n);
        for (int v = 0; v < n; v++) {
            if (v != u) {
                cost_t c = adjacency[u * n + v];
                if (c < INF) {
                    row.push_back({static_cast<uint32_t>(v), c});
                }
            }
        }
        neighbors[u] = std::move(row);
    }

    // ----- HELD-KARP DP -----
    // We skip node0 in subsets, so effectively we have subsets of size up to (n-1).
    const int m = n - 1;
    const size_t dpSize = (static_cast<size_t>(1) << m) * static_cast<size_t>(n);

    // We'll allocate two DP arrays for double-buffering:
    //  dp, dpNext  => cost arrays
    //  parent, parentNext => parent pointers for path reconstruction
    auto dp       = std::unique_ptr<cost_t[]>(new(std::nothrow) cost_t[dpSize]);
    auto dpNext   = std::unique_ptr<cost_t[]>(new(std::nothrow) cost_t[dpSize]);
    auto parent   = std::unique_ptr<int[]>(new(std::nothrow) int[dpSize]);
    auto parentNext = std::unique_ptr<int[]>(new(std::nothrow) int[dpSize]);
    if (!dp || !dpNext || !parent || !parentNext) {
        return TSP_STATUS_OUT_OF_MEMORY;
    }

    // Initialize all DP arrays
    for (size_t i = 0; i < dpSize; i++) {
        dp[i]       = INF;
        dpNext[i]   = INF;
        parent[i]   = -1;
        parentNext[i] = -1;
    }

    // Base cases: single-node subsets {j}, mask = 1<<(j-1) if j>=1
    // dp[mask][j] = cost(0->j)
    for (int j = 1; j < n; j++) {
        cost_t c = adjacency[0 * n + j];
        if (c < INF) {
            size_t mask = 1ULL << (j - 1);
            dp[mask * n + j] = c;
            // parent is -1 meaning came directly from node0
        }
    }

    // Precompute subsetsBySize
    std::vector<std::vector<int>> subsetsBySize(m + 1);
    subsetsBySize.reserve(m + 1);
    {
        const int limit = 1 << m;
        for (int mask = 0; mask < limit; mask++) {
            int c = popcount(mask);
            if (c <= m) {
                subsetsBySize[c].push_back(mask);
            }
        }
    }

    // Main DP: from subsets of size k to subsets of size k+1
    for (int k = 1; k < m; k++) {
        // Reset dpNext, parentNext to “blank” for this layer
        // (because we will fill them from dp)
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (size_t i = 0; i < dpSize; i++) {
            dpNext[i]     = INF;
            parentNext[i] = -1;
        }

        const auto &currentMasks = subsetsBySize[k];

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for (size_t sidx = 0; sidx < currentMasks.size(); sidx++) {
            int mask = currentMasks[sidx];
            // For each j in 'mask' (bit j-1 set => j in subset)
            int subset_j = mask;
            while (subset_j) {
#ifdef __GNUC__
                int jBit = __builtin_ctz(subset_j);
#else
                // fallback:
                // find lowest set bit
                int jBit = 31 - __builtin_clz(subset_j & -subset_j); // or another way
#endif
                subset_j ^= (1 << jBit);
                int j = jBit + 1;  // actual node index

                cost_t cost_j = dp[static_cast<size_t>(mask) * n + j];
                if (cost_j >= INF) {
                    continue;
                }
                // Explore neighbors of j
                for (auto &[nxt, edgeCost] : neighbors[j]) {
                    if (nxt == 0) {
                        // node0 is not visited in the middle
                        continue;
                    }
                    // skip if already visited
                    if (mask & (1 << (nxt - 1))) {
                        continue;
                    }
                    size_t nxtMask = static_cast<size_t>(mask) | (1ULL << (nxt - 1));
                    cost_t newCost = cost_j + edgeCost;
                    size_t idx     = nxtMask * n + nxt;
                    // Relaxation
                    if (newCost < dpNext[idx]) {
                        dpNext[idx] = newCost;
                        parentNext[idx] = j;
                    }
                }
            }
        }

        // Now swap dp <-> dpNext for the next iteration
        dp.swap(dpNext);
        parent.swap(parentNext);
    }

    // Final step: return to node0 from last set
    int fullMask = (1 << m) - 1;
    cost_t best_cost = INF;
    int best_end = -1;
    for (int j = 1; j < n; j++) {
        cost_t c = dp[static_cast<size_t>(fullMask) * n + j];
        if (c >= INF) continue;
        cost_t retC = adjacency[j * n + 0];
        if (retC >= INF) continue;
        cost_t total = c + retC;
        if (total < best_cost) {
            best_cost = total;
            best_end = j;
        }
    }
    if (best_end < 0 || best_cost >= INF) {
        return TSP_STATUS_NO_SOLUTION;
    }

    // Reconstruct
    std::vector<int> pathIndices;
    reconstructPath(parent.get(), n, fullMask, best_end, pathIndices);

    // Final path is 0 -> pathIndices -> 0
    std::vector<int> finalPath;
    finalPath.reserve(pathIndices.size() + 2);
    finalPath.push_back(0);
    for (int node : pathIndices) {
        finalPath.push_back(node);
    }
    finalPath.push_back(0);

    // Convert to original node IDs
    size_t pathLen = finalPath.size();
    auto *sol = new(std::nothrow) uint64_t[pathLen];
    if (!sol) {
        return TSP_STATUS_OUT_OF_MEMORY;
    }
    for (size_t i = 0; i < pathLen; i++) {
        sol[i] = mapper.index_to_id[finalPath[i]];
    }

    output_descriptor->path = sol;
    output_descriptor->num_nodes = pathLen;
    output_descriptor->solution_cost = best_cost;
    return TSP_STATUS_SUCCESS;
}