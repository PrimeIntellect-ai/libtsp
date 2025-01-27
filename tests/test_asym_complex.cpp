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

    std::ifstream inFile("test_data_gen_asym/atsp_test_cases_exact.json");
    ASSERT_TRUE(inFile.is_open()) << "Could not open test_data_gen_asym/atsp_test_cases_exact.json for reading.";

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

        auto edges_array = instance["edges"];
        size_t numEdges = edges_array.size();

        TspInputGraphDescriptor input_desc{};
        input_desc.num_edges = numEdges;
        input_desc.edges = new TspInputGraphEdge[numEdges];

        for (size_t i = 0; i < numEdges; ++i) {
            auto e = edges_array[i];
            input_desc.edges[i].from = e["from"].get<uint64_t>();
            input_desc.edges[i].to = e["to"].get<uint64_t>();
            input_desc.edges[i].cost = e["cost"].get<double>();
        }

        // Solve using tsp_solve
        TspSolutionDescriptor output_desc{};
        TspSolverOptionsDescriptor solver_options{.seed = 0, .num_iterations = 10};

        auto start = std::chrono::high_resolution_clock::now();
        TSPStatus status = tspAsymmetricSolve(&input_desc, &solver_options, &output_desc);
        auto end = std::chrono::high_resolution_clock::now();

        EXPECT_EQ(status, TSP_STATUS_SUCCESS) << "Solver did not return SUCCESS for instance: " << description;

        // Compare with known-good solution
        auto solution_obj = instance["solution"];
        double solution_cost = solution_obj["solution_cost"].get<double>();
        size_t num_nodes = solution_obj["num_nodes"].get<size_t>();

        // assert the numer of nodes in the path is the same as the known-good solution
        EXPECT_EQ(output_desc.num_nodes, num_nodes) << "Path size mismatch for instance: " << description;

        // assert that the path contains the same set of unique node IDs as the known-good solution
        {
            auto knownPath = solution_obj["path"].get<std::vector<uint64_t>>();
            std::unordered_set knownPathSet(knownPath.begin(), knownPath.end());
            std::unordered_set<uint64_t> visited{};
            for (size_t i = 0; i < output_desc.num_nodes; ++i) {
                if (visited.contains(output_desc.tour[i])) {
                    FAIL() << "Path contains duplicate node ID for instance: " << description;
                }
                visited.insert(output_desc.tour[i]);
            }
        }

        // verify if cost is populated correctly for output path
        {
            float pathCost = 0.0f;
            for (size_t i = 0; i < output_desc.num_nodes; ++i) {
                auto from = output_desc.tour[i];
                auto to = output_desc.tour[(i + 1) % output_desc.num_nodes];
                for (size_t j = 0; j < numEdges; ++j) {
                    if (input_desc.edges[j].from == from && input_desc.edges[j].to == to) {
                        pathCost += input_desc.edges[j].cost;
                        break;
                    }
                }
            }
            EXPECT_NEAR(pathCost, output_desc.solution_cost, 1e-6) << "Path cost mismatch for instance: " <<
 description;
        }

        // assert that the cost is within 5% of the known-good solution
        EXPECT_NEAR(output_desc.solution_cost, solution_cost, 0.05 * solution_cost) << "Cost mismatch for instance: " <<
 description;

        std::cout << "Solved \"" << description << "\" with error: " << std::abs(output_desc.solution_cost - solution_cost) << " ("
                << std::fixed << std::setprecision(2) <<
                std::abs(output_desc.solution_cost - solution_cost) / solution_cost * 100 << "%)" << " in "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

        // Cleanup
        delete[] input_desc.edges;
        delete[] output_desc.tour;
    }
}
