#include <gtest/gtest.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdint>
#include <cmath>

#include "json.hpp"

#include <libtsp.h>

/**
 * Checks if two paths are the same up to a cyclic rotation (and ignoring the
 * solver's repetition of the start node at the end).
 */
static bool pathsAreCyclicallyEquivalent(
    const std::vector<uint64_t> &solverPath,
    const std::vector<uint64_t> &knownPath) {
    if (solverPath.size() != knownPath.size()) {
        return false;
    }

    // 1) Try forward cyclic match
    for (size_t start = 0; start < knownPath.size(); ++start) {
        bool match = true;
        for (size_t i = 0; i < knownPath.size(); ++i) {
            if (const size_t idx = (start + i) % knownPath.size(); solverPath[i] != knownPath[idx]) {
                match = false;
                break;
            }
        }
        if (match) {
            return true;
        }
    }

    // 2) Try reversed cyclic match
    // Build the reversed vector of knownPath
    const std::vector reversed(knownPath.rbegin(), knownPath.rend());
    for (size_t start = 0; start < reversed.size(); ++start) {
        bool match = true;
        for (size_t i = 0; i < reversed.size(); ++i) {
            if (const size_t idx = (start + i) % reversed.size(); solverPath[i] != reversed[idx]) {
                match = false;
                break;
            }
        }
        if (match) {
            return true;
        }
    }

    return false;
}

/**
 * Loads the JSON test data from "tsp_data.json" and runs the solver for each TSP instance.
 */
TEST(TspSolverJsonTest, LoadsAndVerifiesExactSolutions) {
    using nlohmann::json;

    std::ifstream inFile("test_data_gen/tsp_test_cases_exact.json");
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
        std::string description = instance["description"].get<std::string>();
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
        TSPStatus status = tsp_solve(&inputDesc, &outputDesc);
        EXPECT_EQ(status, TSP_STATUS_SUCCESS) << "Solver did not return SUCCESS for instance: " << description;

        // Compare with known-good solution
        auto solutionObj = instance["solution"];
        double knownCost = solutionObj["solution_cost"].get<double>();
        size_t knownSize = solutionObj["num_nodes"].get<size_t>();
        auto knownPathJson = solutionObj["path"];
        std::vector<uint64_t> knownPathVec(knownSize);
        for (size_t i = 0; i < knownSize; ++i) {
            knownPathVec[i] = knownPathJson[i].get<uint64_t>();
        }

        // The solver's path presumably includes an extra repeated start node at the end, so let's remove it
        // if it's indeed repeated. (This depends on your solver's implementation!)
        std::vector solverPathVec(outputDesc.path, outputDesc.path + outputDesc.num_nodes);
        if (solverPathVec.size() == knownSize + 1) {
            if (!solverPathVec.empty() && solverPathVec.front() == solverPathVec.back()) {
                solverPathVec.pop_back(); // remove repeated start node
            }
        }

        // Compare path lengths
        EXPECT_EQ(solverPathVec.size(), knownSize);

        // Compare cost with a small epsilon
        EXPECT_NEAR(outputDesc.solution_cost, knownCost, 3e-5);

        // Compare paths (cyclic equivalence)
        EXPECT_TRUE(pathsAreCyclicallyEquivalent(solverPathVec, knownPathVec))
            << "Solver path is not cyclically equivalent to known path for instance: " << description;

        // Cleanup
        delete[] inputDesc.edges;
        delete[] outputDesc.path;
    }
}
