import math
import random
import json
from ortools.linear_solver import pywraplp

def generate_random_points(n, x_range=(0, 100), y_range=(0, 100)):
    """
    Generate n random 2D points within the given ranges.
    Returns a list of (x, y) tuples.
    """
    return [
        (
            random.uniform(x_range[0], x_range[1]),
            random.uniform(y_range[0], y_range[1])
        )
        for _ in range(n)
    ]

def build_distance_matrix(points):
    """
    Given a list of 2D points, compute the pairwise distance matrix,
    where distance_matrix[i][j] = Euclidean distance between points[i] and points[j].
    """
    n = len(points)
    distance_matrix = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                distance_matrix[i][j] = 0.0
            else:
                dx = points[i][0] - points[j][0]
                dy = points[i][1] - points[j][1]
                distance_matrix[i][j] = math.hypot(dx, dy)
    return distance_matrix

def solve_tsp_exact(distance_matrix):
    """
    Solve the TSP exactly using a MIP formulation with OR-Tools.
    Returns (route, cost), where:
      - route is a list of node indices in visit order (starting from 0).
      - cost is the total distance of that route.

    This uses the classic Miller-Tucker-Zemlin (MTZ) constraints to avoid subtours.
    """
    n = len(distance_matrix)
    solver = pywraplp.Solver.CreateSolver("SCIP")  # or "CBC", "BOP", etc.
    if not solver:
        raise RuntimeError("Failed to create MIP solver. Ensure OR-Tools is installed correctly.")

    # Create binary variables x[i,j]: 1 if we go directly from city i to city j, 0 otherwise.
    x = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                x[i, j] = solver.IntVar(0, 1, f"x_{i}_{j}")
            else:
                # Force x[i,i] = 0 with no real variable needed
                # but let's just skip it to keep the dictionary simpler.
                pass

    # Create the "u" variables for MTZ constraints (used to eliminate subtours).
    # We only need them for nodes 1..n-1, but to simplify indexing we'll do 0..n-1.
    # Each u[i] is in [0, n-1], integer.
    u = []
    for i in range(n):
        u.append(solver.IntVar(0, n - 1, f"u_{i}"))

    # Objective: minimize sum of distance_matrix[i][j] * x[i,j]
    objective = []
    for i in range(n):
        for j in range(n):
            if i != j:
                objective.append(distance_matrix[i][j] * x[i, j])
    solver.Minimize(solver.Sum(objective))

    # Constraint 1: Each city i has exactly one outgoing edge.
    for i in range(n):
        solver.Add(solver.Sum([x[i, j] for j in range(n) if j != i]) == 1)

    # Constraint 2: Each city j has exactly one incoming edge.
    for j in range(n):
        solver.Add(solver.Sum([x[i, j] for i in range(n) if i != j]) == 1)

    # MTZ constraint:  u[j] >= u[i] + 1 - n*(1 - x[i,j])  for i != j and i != j != 0
    # Fix u[0] = 0 to define a root (start city).
    solver.Add(u[0] == 0)
    for i in range(1, n):
        for j in range(1, n):
            if i != j:
                # u_j >= u_i + 1 - n*(1 - x[i,j])
                solver.Add(u[j] >= u[i] + 1 - n * (1 - x[i, j]))

    # Solve
    status = solver.Solve()
    if status not in (pywraplp.Solver.OPTIMAL,):
        raise RuntimeError("No feasible TSP solution found within solver constraints or time limit.")

    # Extract the solution route
    # Start at node 0, follow x[0,?] = 1, then x[?,?], etc.
    route = [0]
    total_cost = 0.0
    current_city = 0
    for _ in range(n - 1):
        # Find j where x[current_city,j] = 1
        next_city = None
        for j in range(n):
            if j != current_city and x[current_city, j].solution_value() > 0.5:
                next_city = j
                break
        if next_city is None:
            raise RuntimeError("Unexpected error: No outgoing edge found from city %d" % current_city)
        total_cost += distance_matrix[current_city][next_city]
        route.append(next_city)
        current_city = next_city
    # Complete the cycle back to 0 for total cost
    total_cost += distance_matrix[current_city][0]

    return route, total_cost

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
                    "from": i,  # 64-bit ID if you prefer. For testing, int is fine.
                    "to": j
                })
    return edges_list

def generate_tsp_case(n):
    """
    Generates a TSP case with n nodes:
      1. Creates random points
      2. Builds distance matrix
      3. Solves TSP EXACTLY using MIP + MTZ constraints
      4. Returns a dictionary with 'edges', 'solution.path', 'solution.num_nodes', and 'solution.solution_cost'
    """
    points = generate_random_points(n)
    distance_matrix = build_distance_matrix(points)
    route, route_cost = solve_tsp_exact(distance_matrix)
    edges_list = create_edges_json(distance_matrix)

    return {
        "edges": edges_list,
        "solution": {
            "path": route,           # list of node IDs in visitation order
            "num_nodes": len(route), # should be n
            "solution_cost": route_cost
        }
    }

def main():
    random.seed(123)  # for reproducibility

    sizes = [5, 7, 10, 12, 15, 17, 20, 23, 25, 26, 27, 28, 29, 30]

    all_cases = []
    for size in sizes:
        print(f"Generating TSP case for N={size}...")
        case_data = generate_tsp_case(size)
        print(f"Successfully generated TSP case for N={size}")
        all_cases.append({
            "description": f"Exact TSP instance with N={size}",
            "n": size,
            "edges": case_data["edges"],
            "solution": case_data["solution"]
        })

        # Write out a JSON array of test cases
        with open("tsp_test_cases_exact.json", "w") as f:
            json.dump(all_cases, f, indent=2)

    print("Generated EXACT TSP test cases in 'tsp_test_cases_exact.json'.")

if __name__ == "__main__":
    main()
