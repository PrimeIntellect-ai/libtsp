import json
import math


def main():
    # Load JSON data
    with open("tsp_instance.json", "r") as f:
        data = json.load(f)

    # Basic checks on 'n' and 'solution'
    n = data["n"]
    solution_data = data["solution"]
    path = solution_data["path"]

    # Make sure the solution claims to have n nodes
    assert solution_data["num_nodes"] == n, (
        f"Solution has {solution_data['num_nodes']} nodes but problem states n={n}"
    )

    # Build an adjacency dictionary for cost lookups: edges_dict[(from, to)] = cost
    edges_dict = {}
    for edge in data["edges"]:
        frm, to, cost = edge["from"], edge["to"], edge["cost"]
        edges_dict[(frm, to)] = cost

    # Verify we have exactly n*n (or at least n*(n-1)) edges if it’s a complete TSP
    assert len(edges_dict) == n * (n - 1), "Unexpected number of directed edges"

    # Recompute the solution cost
    recomputed_cost = 0.0
    for i in range(len(path)):
        current_node = path[i]
        # Next node (wrapping around to 0 at the end)
        next_node = path[(i + 1) % len(path)]
        edge_key = (current_node, next_node)

        if edge_key not in edges_dict:
            raise ValueError(f"No edge found for {current_node} -> {next_node}")

        recomputed_cost += edges_dict[edge_key]

    # Assert the recomputed cost matches the provided solution cost within a small epsilon
    stated_cost = solution_data["solution_cost"]
    if not math.isclose(recomputed_cost, stated_cost, rel_tol=1e-9, abs_tol=1e-9):
        raise AssertionError(
            f"Recomputed cost ({recomputed_cost}) does not match "
            f"stated cost ({stated_cost})."
        )

    print("All checks passed. The solution cost is verified to be correct.")


if __name__ == "__main__":
    main()
