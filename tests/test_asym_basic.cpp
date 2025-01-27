#include <gtest/gtest.h>
#include <libtsp.h>
#include <vector>
#include <cstdint>
#include <cmath>

// Helper function to allocate a TspInputGraphDescriptor from a vector of triples (from, to, cost).
static TspInputGraphDescriptor createGraph(const std::vector<std::tuple<uint64_t, uint64_t, double>> &edges) {
    TspInputGraphDescriptor desc{};
    desc.num_edges = edges.size();
    if (!edges.empty()) {
        desc.edges = new TspInputGraphEdge[edges.size()];
        for (size_t i = 0; i < edges.size(); ++i) {
            desc.edges[i].from = std::get<0>(edges[i]);
            desc.edges[i].to = std::get<1>(edges[i]);
            desc.edges[i].cost = std::get<2>(edges[i]);
        }
    }
    return desc;
}

// Helper function to free the allocated arrays within the TspInputGraphDescriptor.
static void freeGraph(TspInputGraphDescriptor &desc) {
    delete[] desc.edges;
    desc.edges = nullptr;
    desc.num_edges = 0;
}

// Helper function to free the path array allocated by tsp_solve.
static void freeOutput(TspSolutionDescriptor &output) {
    delete[] output.tour;
    output.tour = nullptr;
    output.num_nodes = 0;
    output.solution_cost = 0.0;
}

// Each test exercises a different scenario of TSP inputs and verifies
// expected solver behavior.

TEST(TspSolverTest, NullArguments) {
    // With null pointers, we should get TSP_STATUS_ERROR_INVALID_ARG.
    EXPECT_EQ(tspSolveAsymmetric(nullptr, nullptr, nullptr), TSP_STATUS_ERROR_INVALID_ARG);

    TspSolutionDescriptor outputDesc{};
    EXPECT_EQ(tspSolveAsymmetric(nullptr, nullptr, &outputDesc), TSP_STATUS_ERROR_INVALID_ARG);

    constexpr TspSolverOptionsDescriptor optionsDesc{.seed = 42, .num_iterations = 100};
    EXPECT_EQ(tspSolveAsymmetric(nullptr, &optionsDesc, nullptr), TSP_STATUS_ERROR_INVALID_ARG);
}

TEST(TspSolverTest, EmptyEdges) {
    // With no edges, solver should return TSP_STATUS_ERROR_INVALID_ARG.
    TspInputGraphDescriptor graph{};
    graph.edges = nullptr;
    graph.num_edges = 0;

    TspSolutionDescriptor outputDesc{};
    constexpr TspSolverOptionsDescriptor optionsDesc{.seed = 42, .num_iterations = 100};
    EXPECT_EQ(tspSolveAsymmetric(&graph, &optionsDesc, &outputDesc), TSP_STATUS_ERROR_INVALID_ARG);
}

TEST(TspSolverTest, NegativeCostEdge) {
    // Including a negative cost edge should result in TSP_STATUS_ERROR_INVALID_GRAPH.
    const auto inputEdges = std::vector<std::tuple<uint64_t, uint64_t, double>>{
        {1, 2, -1.0}
    };
    TspInputGraphDescriptor graph = createGraph(inputEdges);
    TspSolutionDescriptor outputDesc{};
    constexpr TspSolverOptionsDescriptor optionsDesc{.seed = 42, .num_iterations = 100};
    EXPECT_EQ(tspSolveAsymmetric(&graph, &optionsDesc, &outputDesc), TSP_STATUS_ERROR_INVALID_GRAPH);

    freeOutput(outputDesc);
    freeGraph(graph);
}

TEST(TspSolverTest, TwoNodesNoReturn) {
    // Two nodes with a single direction edge. There's no way to form a tour.
    const auto inputEdges = std::vector<std::tuple<uint64_t, uint64_t, double>>{
        {1, 2, 10.0}
    };
    TspInputGraphDescriptor graph = createGraph(inputEdges);
    TspSolutionDescriptor outputDesc{};
    constexpr TspSolverOptionsDescriptor optionsDesc{.seed = 42, .num_iterations = 100};
    EXPECT_EQ(tspSolveAsymmetric(&graph, &optionsDesc, &outputDesc), TSP_STATUS_ERROR_INVALID_GRAPH);

    freeOutput(outputDesc);
    freeGraph(graph);
}

TEST(TspSolverTest, TwoNodesRoundTrip) {
    // Two nodes with a round trip. The TSP cycle cost is 10 + 10 = 20.
    const auto inputEdges = std::vector<std::tuple<uint64_t, uint64_t, double>>{
        {1, 2, 10.0},
        {2, 1, 10.0}
    };
    TspInputGraphDescriptor graph = createGraph(inputEdges);
    TspSolutionDescriptor outputDesc{};
    constexpr TspSolverOptionsDescriptor optionsDesc{.seed = 42, .num_iterations = 100};
    EXPECT_EQ(tspSolveAsymmetric(&graph, &optionsDesc, &outputDesc), TSP_STATUS_SUCCESS);
    ASSERT_EQ(outputDesc.num_nodes, static_cast<size_t>(2));
    EXPECT_DOUBLE_EQ(outputDesc.solution_cost, 20.0);

    freeOutput(outputDesc);
    freeGraph(graph);
}

TEST(TspSolverTest, DisconnectedThreeNodes) {
    // 3 nodes: 1->2 (5.0), 2->3 (10.0). No return edges to 1. Expect no solution.
    const auto inputEdges = std::vector<std::tuple<uint64_t, uint64_t, double>>{
        {1, 2, 5.0},
        {2, 3, 10.0}
    };
    TspInputGraphDescriptor graph = createGraph(inputEdges);
    TspSolutionDescriptor outputDesc{};
    constexpr TspSolverOptionsDescriptor optionsDesc{.seed = 42, .num_iterations = 100};
    EXPECT_EQ(tspSolveAsymmetric(&graph, &optionsDesc, &outputDesc), TSP_STATUS_ERROR_INVALID_GRAPH);

    freeOutput(outputDesc);
    freeGraph(graph);
}

TEST(TspSolverTest, ThreeNodesFull) {
    // 3 nodes, fully connected with distinct costs.
    // We'll check the path cost is minimal.
    // Suppose 1->2=1, 2->3=2, 3->1=3, 2->1=1, 3->2=2, 1->3=10
    // The best cycle might be 1->2->3->1 with cost=1+2+3=6
    //   or 1->3->2->1 with cost=10+2+1=13
    // So the best is 6.
    const auto inputEdges = std::vector<std::tuple<uint64_t, uint64_t, double>>{
        {1, 2, 1.0},
        {2, 3, 2.0},
        {3, 1, 3.0},

        {2, 1, 1.0},
        {3, 2, 2.0},
        {1, 3, 10.0}
    };
    TspInputGraphDescriptor graph = createGraph(inputEdges);
    TspSolutionDescriptor outputDesc{};
    constexpr TspSolverOptionsDescriptor optionsDesc{.seed = 42, .num_iterations = 100};
    EXPECT_EQ(tspSolveAsymmetric(&graph, &optionsDesc, &outputDesc), TSP_STATUS_SUCCESS);
    EXPECT_NEAR(outputDesc.solution_cost, 6.0, 1e-9);

    freeOutput(outputDesc);
    freeGraph(graph);
}

TEST(TspSolverTest, DuplicateEdges) {
    // Duplicate edges from the same node pair with different costs should result in an error.
    const auto inputEdges = std::vector<std::tuple<uint64_t, uint64_t, double>>{
        {1, 2, 10.0},
        {1, 2, 1.0},
        {2, 1, 1.0},
    };
    TspInputGraphDescriptor graph = createGraph(inputEdges);
    TspSolutionDescriptor outputDesc{};
    constexpr TspSolverOptionsDescriptor optionsDesc{.seed = 42, .num_iterations = 100};
    EXPECT_EQ(tspSolveAsymmetric(&graph, &optionsDesc, &outputDesc), TSP_STATUS_ERROR_INVALID_GRAPH);

    if (outputDesc.tour) {
        freeOutput(outputDesc);
    }
    freeGraph(graph);
}

TEST(TspSolverTest, LargerNoSolution) {
    // 4 nodes in a line, no return edges.
    const auto inputEdges = std::vector<std::tuple<uint64_t, uint64_t, double>>{
        {1, 2, 1.0},
        {2, 3, 1.0},
        {3, 4, 1.0}
        // Missing edges that would allow returning to 1
    };
    TspInputGraphDescriptor graph = createGraph(inputEdges);
    TspSolutionDescriptor outputDesc{};
    constexpr TspSolverOptionsDescriptor optionsDesc{.seed = 42, .num_iterations = 100};
    EXPECT_EQ(tspSolveAsymmetric(&graph, &optionsDesc, &outputDesc), TSP_STATUS_ERROR_INVALID_GRAPH);

    freeOutput(outputDesc);
    freeGraph(graph);
}

TEST(TspSolverTest, FourNodeExample) {
    // Fully connect 4 nodes. We want to ensure the solver returns an actual minimal cycle.
    // We'll connect them in a manner that obviously forms a square: 1->2=1,2->3=1,3->4=1,4->1=1, etc.
    // The minimal cycle could be 1->2->3->4->1 with cost=4.
    // We'll add some extraneous edges with bigger costs, which shouldn't be used.
    const auto inputEdges = std::vector<std::tuple<uint64_t, uint64_t, double>>{
        {1, 2, 1.0}, {2, 3, 1.0}, {3, 4, 1.0}, {4, 1, 1.0}, // cheap cycle
        {2, 1, 10.0}, {3, 2, 10.0}, {4, 3, 10.0}, {1, 4, 10.0}, // expensive return edges

        // Make it fully connected for completeness
        {1, 3, 5.0}, {3, 1, 5.0},
        {2, 4, 5.0}, {4, 2, 5.0}
    };
    TspInputGraphDescriptor graph = createGraph(inputEdges);
    TspSolutionDescriptor outputDesc{};

    constexpr TspSolverOptionsDescriptor optionsDesc{.seed = 42, .num_iterations = 100};
    EXPECT_EQ(tspSolveAsymmetric(&graph, &optionsDesc, &outputDesc), TSP_STATUS_SUCCESS);
    EXPECT_NEAR(outputDesc.solution_cost, 4.0, 1e-9);

    freeOutput(outputDesc);
    freeGraph(graph);
}
