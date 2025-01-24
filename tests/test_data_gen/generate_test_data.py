import math
import random
import json
from ortools.constraint_solver import pywrapcp, routing_enums

def generate_random_points(n, x_range=(0, 100), y_range=(0, 100)):
    """
    Generate n random 2D points in the ranges provided.
    Returns a list of (x, y) tuples.
    """
    return [
        (random.uniform(x_range[0], x_range[1]),
         random.uniform(y_range[0], y_range[1]))
        for _ in range(n)
    ]

def build_distance_matrix(points):
    """
    Given a list of 2D points, compute the pairwise distance matrix.
    distance_matrix[i][j] = Euclidean distance between points[i] and points[j].
    """
    n = len(points)
    distance_matrix = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                distance_matrix[i][j] = 0.0
            else:
                distance_matrix[i][j] = math.dist(points[i], points[j])
    return distance_matrix

def solve_tsp_ortools(distance_matrix):
    """
    Solve the TSP using OR-Tools and return the best route + minimal cost.

    Returns:
      - A list of node indices representing the route (excluding the repeated start at the end).
      - The total cost of that route.
    """
    n = len(distance_matrix)
    manager = pywrapcp.RoutingIndexManager(n, 1, 0)  # 1 vehicle, start node = 0
    routing = pywrapcp.RoutingModel(manager)

    def distance_callback(from_index, to_index):
        # Convert routing variable Index to distance matrix node index.
        from_node = manager.IndexToNode(from_index)
        to_node = manager.IndexToNode(to_index)
        return int(distance_matrix[from_node][to_node])  # OR-Tools demands int or long

    transit_callback_index = routing.RegisterTransitCallback(distance_callback)
    routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)

    # We want to find a route that visits all nodes exactly once.
    # Set search parameters
    search_parameters = pywrapcp.DefaultRoutingSearchParameters()
    search_parameters.first_solution_strategy = (
        routing_enums.FirstSolutionStrategy.PATH_CHEAPEST_ARC
    )
    # Use a more advanced strategy (or extended search) to get optimal solutions.
    search_parameters.local_search_metaheuristic = (
        routing_enums.LocalSearchMetaheuristic.GUIDED_LOCAL_SEARCH
    )
    search_parameters.time_limit.seconds = 30  # Increase if needed
    search_parameters.log_search = False

    solution = routing.SolveWithParameters(search_parameters)
    if solution:
        # Extract the route
        index = routing.Start(0)
        route = []
        route_cost = 0
        while not routing.IsEnd(index):
            route.append(manager.IndexToNode(index))
            next_index = solution.Value(routing.NextVar(index))
            route_cost += distance_matrix[manager.IndexToNode(index)][manager.IndexToNode(next_index)]
            index = next_index

        # route_cost is a float. The route is the order in which nodes are visited (starting at 0).
        return route, route_cost
    else:
        raise RuntimeError("No solution found by OR-Tools.")

def create_edges_json(distance_matrix):
    """
    Create the JSON array of edges (with fields: cost, from, to).
    Each pair of distinct nodes i, j is stored as an edge from i to j.
    """
    edges_list = []
    n = len(distance_matrix)
    for i in range(n):
        for j in range(n):
            if i != j:
                edges_list.append({
                    "cost": distance_matrix[i][j],
                    "from": i,  # Or use a 64-bit ID scheme if you prefer
                    "to": j
                })
    return edges_list


def generate_tsp_case(n):
    """
    Generates a TSP case for n nodes:
    - Creates random points
    - Builds distance matrix
    - Solves TSP
    - Returns a dictionary with 'edges', 'path', 'num_nodes', and 'solution_cost'
    """
    points = generate_random_points(n)
    distance_matrix = build_distance_matrix(points)

    route, route_cost = solve_tsp_ortools(distance_matrix)
    edges_list = create_edges_json(distance_matrix)

    # The solution path is an array of node IDs (64-bit ints are fine; here we just use int).
    # The cost is a non-negative double.
    tsp_case = {
        "edges": edges_list,
        "solution": {
            "path": route,
            "num_nodes": len(route),
            "solution_cost": route_cost
        }
    }
    return tsp_case


def main():
    # Example: generate multiple TSP test cases of various sizes
    random.seed(123)  # for reproducibility

    # You can adjust these sizes or just do a single N
    sizes = [5, 7, 10]  # up to 25 if you like
    all_cases = []

    for size in sizes:
        case_data = generate_tsp_case(size)
        all_cases.append({
            "description": f"TSP instance with N={size}",
            "n": size,
            "edges": case_data["edges"],
            "solution": case_data["solution"]
        })

    # Write out a JSON array of test cases
    with open("tsp_test_cases.json", "w") as f:
        json.dump(all_cases, f, indent=2)

    print("Generated TSP test cases in 'tsp_test_cases.json'.")

if __name__ == "__main__":
    main()
