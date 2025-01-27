#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>
#include <libtsp.h>

static TspInputGraphDescriptor createRandomFullyConnectedGraph(const size_t n, const double max_cost) {
    TspInputGraphDescriptor desc{};
    if (n == 0) {
        return desc;
    }

    std::mt19937_64 rng{std::random_device{}()};
    std::uniform_real_distribution dist(1.0, max_cost);

    const size_t num_edges = n * (n - 1);
    desc.edges = new TspInputGraphEdge[num_edges];
    desc.num_edges = num_edges;

    size_t index = 0;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            if (i == j) {
                continue;
            }
            desc.edges[index].from = i;
            desc.edges[index].to   = j;
            desc.edges[index].cost = dist(rng);
            index++;
        }
    }
    return desc;
}

static void freeGraph(TspInputGraphDescriptor& desc) {
    delete[] desc.edges;
    desc.edges = nullptr;
    desc.num_edges = 0;
}

static void freeOutput(TspSolutionDescriptor& output) {
    delete[] output.tour;
    output.tour = nullptr;
    output.num_nodes = 0;
    output.solution_cost = 0.0;
}

int main() {
    for (const std::vector<size_t> test_sizes{5, 7, 9, 11, 13, 15, 17, 20, 22, 24, 25}; const auto n : test_sizes) {
        TspInputGraphDescriptor graph = createRandomFullyConnectedGraph(n, 100.0);

        constexpr size_t iterations = 3;
        long long total_ns = 0; 

        for (size_t i = 0; i < iterations; i++) {
            TspSolutionDescriptor outputDesc{};

            auto start = std::chrono::steady_clock::now();
            if (const TSPStatus status = tspSolveSymmetric(&graph, {}, &outputDesc); status != TSP_STATUS_SUCCESS) {
                std::cerr << "Error: tsp_solve returned status " << status << std::endl;
                freeOutput(outputDesc);
                freeGraph(graph);
                return 1;
            }
            auto end = std::chrono::steady_clock::now();

            const long long duration_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
            total_ns += duration_ns;

            freeOutput(outputDesc);
        }

        freeGraph(graph);

        const double avg_ms = (static_cast<double>(total_ns) / iterations) / 1e6;

        // Print a summary for this n
        std::cout << "N = " << n
                  << ", Average Solve Time = " << std::fixed << std::setprecision(3) << avg_ms
                  << " ms (over " << iterations << " runs)" << std::endl;
    }

    return 0;
}
