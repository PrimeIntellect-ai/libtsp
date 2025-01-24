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

#include <libtsp.h> // Adjust if needed

static constexpr cost_t INF = std::numeric_limits<cost_t>::infinity();

namespace {
    // Maps arbitrary node IDs to [0..n-1].
    struct NodeMapper {
        std::vector<uint64_t> index_to_id;
        std::map<uint64_t, int> id_to_index;

        void build(const TspInputGraphDescriptor *graph) {
            index_to_id.reserve(graph->num_edges * 2ULL);
            for (size_t i = 0; i < graph->num_edges; ++i) {
                auto from_id = graph->edges[i].from;
                auto to_id = graph->edges[i].to;
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

    // A fast popcount for 32-bit. For bigger sets up to n=25, still fits in 32 bits.
    int popcount(unsigned x) {
#ifdef __GNUC__
        return __builtin_popcount(x);
#else
    // fallback if non-GCC:
    int c=0;
    while (x) {
        x &= (x-1);
        c++;
    }
    return c;
#endif
    }

    // Reconstruct path from DP parent info.
    // dp[mask*n + j] = cost of traveling from node0 -> ... -> j with visited set = mask.
    // parent[mask*n + j] = previous node in the best path.
    void reconstructPath(
        const int *parent,
        const int n,
        const int fullMask,
        const int endNode,
        std::vector<int> &pathOut) {
        pathOut.clear();
        int curMask = fullMask;
        int curNode = endNode;
        while (true) {
            pathOut.push_back(curNode);
            const int p = parent[static_cast<size_t>(curMask) * n + curNode];
            if (p < 0) {
                // means we came from node 0
                break;
            }
            // remove curNode from mask
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
                               TspOutputGraphDescriptor *output_descriptor) {
    // Basic checks
    if (!graph || !output_descriptor) {
        return TSP_STATUS_ERROR_INVALID_ARG;
    }
    if (!graph->edges || graph->num_edges == 0) {
        return TSP_STATUS_ERROR_INVALID_ARG;
    }

    NodeMapper mapper;
    mapper.build(graph);
    int n = static_cast<int>(mapper.size());
    if (n == 0) {
        return TSP_STATUS_NO_SOLUTION;
    }
    if (n == 1) {
        // Only one node => cost=0 path.
        output_descriptor->num_nodes = 1;
        output_descriptor->path = new(std::nothrow) uint64_t[1];
        if (!output_descriptor->path) {
            return TSP_STATUS_OUT_OF_MEMORY;
        }
        output_descriptor->path[0] = mapper.index_to_id[0];
        output_descriptor->solution_cost = 0.0;
        return TSP_STATUS_SUCCESS;
    }

    // Build full adjacency matrix, adjacency[u*n + v] = cost(u->v) or INF if missing.
    std::vector adjacency(n * n, INF);
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

    // Quick handle n=2
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

    // Edge-Dominance Elimination
    // For each edge u->v, check if cost(u->v) > cost(u->w) + cost(w->v) for some w != u,v.
    // If yes, eliminate (u->v) by setting adjacency[u*n + v] = INF.
    // This is O(n^3). For n <= 25, it's borderline but often still acceptable.

    // We'll do: for (u in [0..n-1]) for (v in [0..n-1]) if cost(u->v)<INF then check for any w.
    for (uint32_t u = 0; u < n; u++) {
        for (uint32_t v = 0; v < n; v++) {
            cost_t uv = adjacency[u * n + v];
            if (u != v && uv < INF) {
                // check dominance
                for (int w = 0; w < n; w++) {
                    if (w == u || w == v) continue;
                    cost_t uw = adjacency[u * n + w];
                    cost_t wv = adjacency[w * n + v];
                    if (uw < INF && wv < INF) {
                        if (uv > uw + wv) {
                            adjacency[u * n + v] = INF;
                            break; // no need to keep checking once eliminated
                        }
                    }
                }
            }
        }
    }

    // Build Sparse Adjacency Lists from the pruned adjacency matrix
    // neighbors[u] = vector of (v, cost(u->v)) for edges not INF.
    std::vector<std::vector<IntCostPair>> neighbors(n);
    neighbors.shrink_to_fit();
    for (uint32_t u = 0; u < n; u++) {
        std::vector<IntCostPair> row;
        row.reserve(n);
        for (uint32_t v = 0; v < n; v++) {
            if (cost_t c = adjacency[u * n + v]; v != u && c < INF) {
                row.push_back({v, c});
            }
        }
        neighbors[u] = std::move(row);
    }

    // Held-Karp DP, skipping node0 from subset
    int m = n - 1;
    size_t dpSize = (static_cast<size_t>(1) << m) * static_cast<size_t>(n);
    auto dp = std::unique_ptr<cost_t[]>(new(std::nothrow) cost_t[dpSize]);
    auto parent = std::unique_ptr<int[]>(new(std::nothrow) int[dpSize]);
    if (!dp || !parent) {
        return TSP_STATUS_OUT_OF_MEMORY;
    }
    for (size_t i = 0; i < dpSize; i++) {
        dp[i] = INF;
        parent[i] = -1;
    }

    // Base cases: single-node subsets {j}, mask = (1 << (j-1)) if j>=1
    // dp[mask][j] = cost(0->j)
    for (int j = 1; j < n; j++) {
        // cost from 0->j
        // we'll find it in adjacency[0*n + j], or use the neighbors[0] list
        cost_t c = adjacency[0 * n + j];
        if (c < INF) {
            size_t mask = 1ULL << (j - 1);
            dp[mask * n + j] = c;
        }
    }

    // We'll process subsets by cardinality
    // Precompute subsets by size to ensure we do not overwrite data for the same subset size.
    std::vector<std::vector<int>> subsetsBySize(m + 1);
    subsetsBySize.reserve(m + 1); {
        int limit = 1 << m;
        for (int mask = 0; mask < limit; mask++) {
            if (int c = popcount(mask); c <= m) {
                subsetsBySize[c].push_back(mask);
            }
        }
    }

    // DP main loop
    for (int k = 1; k < m; k++) {
        // For each subset of size k
        for (const auto &levelMasks = subsetsBySize[k]; int mask: levelMasks) {
            // For each 'end' node j in subset => bit j-1 in mask
            int subset_j = mask;
            while (subset_j) {
                int jBit = __builtin_ctz(subset_j);
                subset_j ^= (1 << jBit);
                int j = jBit + 1; // node index
                cost_t cost_j = dp[static_cast<size_t>(mask) * n + j];
                if (cost_j >= INF) continue;

                // Now iterate only over the neighbors of j
                // If neighbor = kNode, skip if it's in mask => already visited
                for (auto &[fst, snd]: neighbors[j]) {
                    int nxt = fst;
                    cost_t edgeCost = snd;
                    // node nxt is in [0..n-1], skip if nxt=0 or if nxt is in subset
                    if (nxt == 0) {
                        continue; // can't visit node0 in the middle
                    }
                    // check if bit (nxt-1) in mask is set
                    if (mask & (1 << (nxt - 1))) {
                        // already visited
                        continue;
                    }
                    size_t nxtMask = static_cast<size_t>(mask) | (1ULL << (nxt - 1));
                    cost_t newCost = cost_j + edgeCost;
                    if (cost_t &dpRef = dp[nxtMask * n + nxt]; newCost < dpRef) {
                        dpRef = newCost;
                        parent[nxtMask * n + nxt] = j;
                    }
                }
            }
        }
    }

    // Final step: return to node0
    int fullMask = (1 << m) - 1;
    cost_t best_cost = INF;
    int best_end = -1;
    for (int j = 1; j < n; j++) {
        cost_t c = dp[static_cast<size_t>(fullMask) * n + j];
        if (c >= INF) continue;
        // cost j->0
        cost_t retC = adjacency[j * n + 0];
        if (retC >= INF) continue;
        if (cost_t total = c + retC; total < best_cost) {
            best_cost = total;
            best_end = j;
        }
    }

    if (best_end < 0 || best_cost >= INF) {
        return TSP_STATUS_NO_SOLUTION;
    }

    // Reconstruct path
    std::vector<int> pathIndices;
    reconstructPath(parent.get(), n, fullMask, best_end, pathIndices);

    // pathIndices is [1..some..end], so final path: 0 -> pathIndices -> 0
    // We'll create finalPath with those
    std::vector<int> finalPath;
    finalPath.reserve(pathIndices.size() + 2);
    finalPath.push_back(0);
    for (int node: pathIndices) {
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
