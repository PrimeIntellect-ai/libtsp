#include <gtest/gtest.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdint>
#include <cmath>

#include "json.hpp"

#include <libtsp.h>
#include <unordered_set>

/**
 * Loads the JSON test data from "tsp_data.json" and runs the solver for each TSP instance.
 * Asserts that the solver return cost within 30% of the known-good solution.
 */
TEST(TspSolverJsonTest, CompareCosts) {
    using nlohmann::json;

    std::ifstream inFile("test_data_gen_sym/tsp_test_cases_exact.json");
    ASSERT_TRUE(inFile.is_open()) << "Could not open test_data_gen/tsp_test_cases_exact.json for reading.";

    json testData;
    inFile >> testData;
    inFile.close();

    // For each test instance, build the TspInputGraphDescriptor, run the solver, and compare results
    for (auto &instance: testData) {
        // The file has the following structure per instance:
        // {
        //   "description": "Exact TSP instance with N=5",
        //   "n": 5,
        //   "edges": [...],
        //   "solution": {
        //     "path": [...],
        //     "num_nodes": 5,
        //     "solution_cost": 188.762845074138
        //   }
        // }
        auto description = instance["description"].get<std::string>();
        SCOPED_TRACE("Testing TSP instance: " + description);

        auto edgesArray = instance["edges"];
        size_t numEdges = edgesArray.size();

        TspInputGraphDescriptor inputDesc{};
        inputDesc.num_edges = numEdges;
        inputDesc.edges = new TspInputGraphEdge[numEdges];

        for (size_t i = 0; i < numEdges; ++i) {
            auto e = edgesArray[i];
            inputDesc.edges[i].from = e["from"].get<uint64_t>();
            inputDesc.edges[i].to = e["to"].get<uint64_t>();
            inputDesc.edges[i].cost = e["cost"].get<double>();
        }

        // Solve using tsp_solve
        TspOutputGraphDescriptor outputDesc{};
        TSPStatus status = tspSolveSymmetric(&inputDesc, {}, &outputDesc);
        EXPECT_EQ(status, TSP_STATUS_SUCCESS) << "Solver did not return SUCCESS for instance: " << description;

        // Compare with known-good solution
        auto solutionObj = instance["solution"];
        double knownCost = solutionObj["solution_cost"].get<double>();
        size_t knownSize = solutionObj["num_nodes"].get<size_t>();

        // assert the numer of nodes in the path is the same as the known-good solution
        EXPECT_EQ(outputDesc.num_nodes, knownSize) << "Path size mismatch for instance: " << description;

        // assert that the path contains the same set of unique node IDs as the known-good solution
        {
            auto knownPath = solutionObj["path"].get<std::vector<uint64_t>>();
            std::unordered_set knownPathSet(knownPath.begin(), knownPath.end());
            std::unordered_set<uint64_t> visited{};
            for (size_t i = 0; i < outputDesc.num_nodes; ++i) {
                if (visited.contains(outputDesc.tour[i])) {
                    FAIL() << "Path contains duplicate node ID for instance: " << description;
                }
                visited.insert(outputDesc.tour[i]);
            }
        }

        // verify if cost is populated correctly for output path
        {
            double pathCost = 0.0;
            for (size_t i = 0; i < outputDesc.num_nodes; ++i) {
                auto from = outputDesc.tour[i];
                auto to = outputDesc.tour[(i + 1) % outputDesc.num_nodes];
                for (size_t j = 0; j < numEdges; ++j) {
                    if (inputDesc.edges[j].from == from && inputDesc.edges[j].to == to) {
                        pathCost += inputDesc.edges[j].cost;
                        break;
                    }
                }
            }
            EXPECT_NEAR(pathCost, outputDesc.solution_cost, 1e-9) << "Path cost mismatch for instance: " << description;
        }


        // assert that the cost is within 15% of the known-good solution
        EXPECT_NEAR(outputDesc.solution_cost, knownCost, 0.15 * knownCost) << "Cost mismatch for instance: " << description;

        // Cleanup
        delete[] inputDesc.edges;
        delete[] outputDesc.tour;
    }
}
