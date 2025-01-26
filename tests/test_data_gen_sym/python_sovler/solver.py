import json
import math
import random
import time
from enum import IntEnum
from dataclasses import dataclass, field
from typing import List, Tuple, Optional

import numpy as np


# =========================================================
# TSPStatus (mirroring the TSPResult enum)
# =========================================================
class TSPStatus(IntEnum):
    """
    Mirrors the enum TSPResult in the provided C code.
    """
    TSP_STATUS_SUCCESS = 0  # Operation completed successfully
    TSP_STATUS_ERROR_INVALID_GRAPH = 1
    TSP_STATUS_ERROR_INVALID_ARG = 2
    TSP_STATUS_NO_SOLUTION = 3
    TSP_STATUS_OUT_OF_MEMORY = 4
    TSP_STATUS_ERROR_INTERNAL = 5


# =========================================================
# TspInputGraphEdge structure
# =========================================================
@dataclass
class TspInputGraphEdge:
    """
    Mirrors the struct TspInputGraphEdge in the provided C code.
    """
    cost: float = 0.0
    from_id: int = 0
    to_id: int = 0


# =========================================================
# TspInputGraphDescriptor structure
# =========================================================
@dataclass
class TspInputGraphDescriptor:
    """
    Mirrors the struct TspInputGraphDescriptor in the provided C code.
    """
    edges: List[TspInputGraphEdge] = field(default_factory=list)
    num_edges: int = 0


# =========================================================
# TspOutputGraphDescriptor structure
# =========================================================
@dataclass
class TspOutputGraphDescriptor:
    """
    Mirrors the struct TspOutputGraphDescriptor in the provided C code.
    """
    path: List[int] = field(default_factory=list)
    num_nodes: int = 0
    solution_cost: float = 0.0


# ========== The TSP Solver ==========
def tsp_solve(
        graph: TspInputGraphDescriptor,
        output_descriptor: TspOutputGraphDescriptor
) -> TSPStatus:
    """
    A TSP solver that:
      1) Builds an adjacency matrix from the input descriptor.
      2) Constructs an initial tour using a Nearest-Neighbor approach.
      3) Improves that tour with a 2-Opt local search pass.
      4) Tries a random double-bridge (4-opt) move to escape local minima,
         then re-runs 2-Opt if it improved.

    This is still simpler than the full Lin-Kernighan, but
    demonstrates an extra "double-bridge" step from that family.
    """

    # ------------------------------------------------------------------
    # 1) Basic Input & Graph Checks
    # ------------------------------------------------------------------
    if not graph or not isinstance(graph, TspInputGraphDescriptor):
        return TSPStatus.TSP_STATUS_ERROR_INVALID_ARG
    if not graph.edges or graph.num_edges < 1:
        return TSPStatus.TSP_STATUS_ERROR_INVALID_ARG

    # Identify the number of nodes (assuming node IDs in [0..N-1])
    max_node_id = -1
    for e in graph.edges:
        max_node_id = max(max_node_id, e.from_id, e.to_id)
    num_nodes = max_node_id + 1
    if num_nodes < 2:
        return TSPStatus.TSP_STATUS_ERROR_INVALID_GRAPH

    # Build adjacency matrix; default to infinity
    adjacency = [[math.inf] * num_nodes for _ in range(num_nodes)]
    for e in graph.edges:
        adjacency[e.from_id][e.to_id] = e.cost

    # Basic connectivity check
    for i in range(num_nodes):
        if all(math.isinf(adjacency[i][j]) for j in range(num_nodes) if j != i):
            return TSPStatus.TSP_STATUS_ERROR_INVALID_GRAPH

    # ------------------------------------------------------------------
    # 2) Helper Functions: Compute Tour Cost, 2-Opt, etc.
    # ------------------------------------------------------------------
    def compute_tour_cost(p: List[int]) -> float:
        # Sum edges p[i] -> p[i+1], plus the cycle edge p[-1] -> p[0]
        total = 0.0
        for i in range(len(p) - 1):
            total += adjacency[p[i]][p[i+1]]
        total += adjacency[p[-1]][p[0]]
        return total

    def two_opt_swap(p: List[int], i: int, j: int) -> None:
        """
        Reverses the sub-path p[i+1..j] in-place. Standard 2-opt edge swap.
        """
        p[i+1:j+1] = reversed(p[i+1:j+1])

    def two_opt_improve(p: List[int]) -> List[int]:
        """
        Runs repeated first-improvement 2-Opt until no improvement is found.
        Modifies the path in-place, then returns it.
        """
        n = len(p)
        improved = True
        while improved:
            improved = False
            for i in range(n - 1):
                for j in range(i + 2, n):
                    # Avoid a swap that cuts the 'wrap-around' edge if j == n-1 and i == 0
                    if j == n - 1 and i == 0:
                        continue
                    # Current edges: (p[i], p[i+1]) and (p[j], p[(j+1)%n])
                    old_cost = (adjacency[p[i]][p[i+1]] +
                                adjacency[p[j]][p[(j+1)%n]])
                    new_cost = (adjacency[p[i]][p[j]] +
                                adjacency[p[i+1]][p[(j+1)%n]])
                    if new_cost < old_cost:
                        two_opt_swap(p, i, j)
                        improved = True
                        break
                if improved:
                    break
        return p

    def double_bridge_move(p: List[int], i: int, j: int, k: int, l: int) -> List[int]:
        """
        Executes a 'double-bridge' 4-opt style move on the cycle path.
        The path is treated as a ring, but for simplicity, we treat it as a
        linear array ignoring the final link p[-1]->p[0] while splicing.

        The standard double-bridge rearrangement can be described as splitting
        the cycle into five segments:
          A = p[0..i]
          B = p[i+1..j]
          C = p[j+1..k]
          D = p[k+1..l]
          E = p[l+1..end]
        We then reassemble them as:
          newPath = A + D + C + B + E

        This typically yields a different cycle that might jump out
        of 2-Opt local minima.
        """
        # Slice carefully. Python slicing excludes the end index
        A = p[:i+1]
        B = p[i+1:j+1]
        C = p[j+1:k+1]
        D = p[k+1:l+1]
        E = p[l+1:]
        return A + D + C + B + E

    def try_random_double_bridge(p: List[int], tries=50) -> bool:
        """
        Attempt up to 'tries' random double-bridge transformations:
        - Randomly pick 4 distinct indices i<j<k<l
        - If the new cycle is better, accept it & return True
        - If no improvement found, return False

        For demonstration, we only do a single acceptance: once we find
        an improved tour, we keep it and stop.
        """
        n = len(p)
        if n < 5:
            return False  # Not enough nodes for 4 distinct breakpoints

        best_cost = compute_tour_cost(p)
        original_path = p[:]

        for _ in range(tries):
            # Randomly pick i<j<k<l
            i = random.randint(0, n-4)
            j = random.randint(i+1, n-3)
            k = random.randint(j+1, n-2)
            l = random.randint(k+1, n-1)
            new_path = double_bridge_move(p, i, j, k, l)
            new_cost = compute_tour_cost(new_path)
            if new_cost < best_cost:
                # Found an improvement
                p[:] = new_path  # copy back to 'p' in-place
                return True
        # No improvement found
        p[:] = original_path  # restore
        return False

    # ------------------------------------------------------------------
    # 3) Build an Initial Tour: Nearest Neighbor
    # ------------------------------------------------------------------
    start_node = 0
    unvisited = set(range(num_nodes))
    unvisited.remove(start_node)
    path = [start_node]
    curr_node = start_node

    while unvisited:
        best_next = None
        best_cost = math.inf
        for nxt in unvisited:
            cost = adjacency[curr_node][nxt]
            if cost < best_cost:
                best_cost = cost
                best_next = nxt
        if best_next is None or math.isinf(best_cost):
            return TSPStatus.TSP_STATUS_NO_SOLUTION
        path.append(best_next)
        unvisited.remove(best_next)
        curr_node = best_next

    # ------------------------------------------------------------------
    # 4) Improve the Tour: 2-Opt
    # ------------------------------------------------------------------
    two_opt_improve(path)

    # ------------------------------------------------------------------
    # 5) Extra Refinement: Try Double-Bridge Moves
    # ------------------------------------------------------------------
    improved = try_random_double_bridge(path, tries=50)
    if improved:
        # If we found a beneficial double-bridge, do 2-Opt again
        two_opt_improve(path)

    # ------------------------------------------------------------------
    # 6) Finalize Output
    # ------------------------------------------------------------------
    final_cost = compute_tour_cost(path)
    output_descriptor.num_nodes = num_nodes
    output_descriptor.path = path
    output_descriptor.solution_cost = final_cost
    return TSPStatus.TSP_STATUS_SUCCESS



def parse_tsplib_with_solution(
        tsp_file_path: str,
        sol_file_path: Optional[str] = None
) -> Tuple[TspInputGraphDescriptor, List[int], float]:
    """
    Parses a TSPLIB .tsp file into a TspInputGraphDescriptor using vectorized NumPy operations
    for fast adjacency matrix construction. If a .sol file is provided, we also parse
    the solution path and compute its total cost.

    :param tsp_file_path: Path to the .tsp file in TSPLIB format.
    :param sol_file_path: Optional path to a .sol file containing a known solution path.
    :return: A tuple of (graph_descriptor, solution_path, solution_cost):
       - graph_descriptor: TspInputGraphDescriptor with edges for a full graph.
       - solution_path: List of node indices (from the .sol file) or empty if none.
       - solution_cost: The total cost of that path (or 0.0 if no .sol was provided).
    """
    # ---------------------------
    # 1) Parse the .tsp file
    # ---------------------------
    dimension = None
    edge_weight_type = None

    # We'll store node coordinates in a NumPy array.
    # But we can’t allocate until we know 'dimension'.
    node_coords_list = []

    with open(tsp_file_path, 'r') as f:
        lines = f.readlines()

    idx = 0
    while idx < len(lines):
        line = lines[idx].strip()
        idx += 1

        if not line or line.upper().startswith("COMMENT"):
            continue

        if line.upper().startswith("NAME:"):
            # e.g. "NAME: berlin52"
            pass  # We don't necessarily need 'name' here, so we can skip storing it.
        elif line.upper().startswith("TYPE:"):
            # e.g. "TYPE: TSP"
            pass
        elif line.upper().startswith("DIMENSION:"):
            dimension = int(line.split(":", 1)[1].strip())
        elif line.upper().startswith("EDGE_WEIGHT_TYPE:"):
            edge_weight_type = line.split(":", 1)[1].strip().upper()
        elif line.upper().startswith("NODE_COORD_SECTION"):
            if dimension is None:
                raise ValueError("DIMENSION not found before NODE_COORD_SECTION.")

            for _ in range(dimension):
                coord_line = lines[idx].strip()
                idx += 1
                parts = coord_line.split()
                # parts[0] = node_id, parts[1] = x, parts[2] = y
                if len(parts) < 3:
                    raise ValueError(f"Invalid NODE_COORD line: {coord_line}")
                x, y = float(parts[1]), float(parts[2])
                node_coords_list.append((x, y))
        elif line.upper().startswith("EOF"):
            break

    if dimension is None or len(node_coords_list) != dimension:
        raise ValueError("TSPLIB file parsing error: dimension or node_coords are missing/inconsistent.")

    if edge_weight_type is None:
        edge_weight_type = "EUC_2D"  # default assumption if not specified

    # ---------------------------
    # 2) Build the adjacency matrix with NumPy
    # ---------------------------
    # Convert the collected list of coords to a NumPy array of shape (N, 2).
    node_coords = np.array(node_coords_list, dtype=float)  # shape: (dimension, 2)

    # We can compute all pairwise distances using broadcasting:
    # delta[i, j, :] = node_coords[i, :] - node_coords[j, :]
    # then dist[i, j] = sqrt( delta[i, j, 0]^2 + delta[i, j, 1]^2 )
    # shape of delta: (dimension, dimension, 2)
    delta = node_coords[:, np.newaxis, :] - node_coords[np.newaxis, :, :]
    distances = np.sqrt(np.sum(delta ** 2, axis=2))  # shape (dimension, dimension)

    if edge_weight_type == "EUC_2D":
        # Round according to TSPLIB spec: floor(dist + 0.5)
        # Convert to float in case dimension is large.
        distances = np.floor(distances + 0.5)

    # Set the diagonal to infinity (no self-loops in TSP).
    np.fill_diagonal(distances, np.inf)

    # Build TspInputGraphEdge list from the matrix
    # We'll only add edges for i != j
    edges = []
    for i in range(dimension):
        for j in range(dimension):
            if i != j:
                cost = float(distances[i, j])  # ensure it's a Python float
                edges.append(TspInputGraphEdge(cost=cost, from_id=i, to_id=j))

    graph_descriptor = TspInputGraphDescriptor(edges=edges, num_edges=len(edges))

    # ---------------------------
    # 3) Parse the .sol file (optional) and compute cost
    # ---------------------------
    solution_path: List[int] = []
    solution_cost: float = 0.0

    if sol_file_path:
        with open(sol_file_path, 'r') as f_sol:
            sol_lines = f_sol.read().strip().split()

        # The first integer in the .sol file is typically the dimension
        sol_dimension = int(sol_lines[0])
        if sol_dimension != dimension:
            raise ValueError(
                f"Mismatch in dimension between .tsp ({dimension}) and .sol ({sol_dimension})."
            )

        # The remaining lines/numbers are the path
        raw_path = sol_lines[1:]
        if len(raw_path) != dimension:
            raise ValueError("Solution path in .sol does not match the dimension.")

        solution_path = [int(x) for x in raw_path]

        # Compute the total cost of traveling the path in a cycle
        route_cost = 0.0
        for i in range(dimension - 1):
            route_cost += distances[solution_path[i], solution_path[i + 1]]
        # Add cost for closing the loop
        route_cost += distances[solution_path[-1], solution_path[0]]

        solution_cost = route_cost

    return graph_descriptor, solution_path, solution_cost


def parse_tsp_json(json_file_path: str) -> List[Tuple[TspInputGraphDescriptor, List[int], float]]:
    with open(json_file_path, 'r') as f:
        data = json.load(f)

    results = []
    for instance in data:
        # Parse edges
        edges_list = instance.get("edges", [])
        edges = [
            TspInputGraphEdge(
                cost=edge["cost"],
                from_id=edge["from"],
                to_id=edge["to"]
            )
            for edge in edges_list
        ]

        # Create TspInputGraphDescriptor
        graph_descriptor = TspInputGraphDescriptor(
            edges=edges,
            num_edges=len(edges)
        )

        # Parse the solution information
        solution_info = instance.get("solution", {})
        solution_path = solution_info.get("path", [])
        solution_cost = solution_info.get("solution_cost", 0.0)

        # Collect them in a tuple
        results.append((graph_descriptor, solution_path, solution_cost))

    return results


def run_test_case(idx: int, graph: TspInputGraphDescriptor, expected_path: List[int], expected_cost: float):
    output_descriptor = TspOutputGraphDescriptor()
    start_time = time.time()
    result = tsp_solve(graph, output_descriptor)
    assert result == TSPStatus.TSP_STATUS_SUCCESS, f"Test case {idx + 1} failed with status {result}"
    end_time = time.time()

    print(f"Test Case {idx + 1}:")
    print(f"  Number of edges: {graph.num_edges}")
    print(f"  Number of nodes: {output_descriptor.num_nodes}")
    print(f"  Expected Path: {expected_path}")
    print(f"  Output Path: {output_descriptor.path}")
    print(f"  Expected Cost: {expected_cost}")
    print(f"  Output Cost: {output_descriptor.solution_cost}")
    print(f"  Time taken: {end_time - start_time:.6f} seconds")
    print()


# =========================================================
# Main function to run the TSP solver on test cases
# =========================================================
def main_1():
    tsp_file = "tsp_cases.json"
    test_cases = parse_tsp_json(tsp_file)

    for idx, (graph, expected_path, expected_cost) in enumerate(test_cases):
        run_test_case(idx, graph, expected_path, expected_cost)


if __name__ == "__main__":
    main_1()
    pass


def main_2():
    berlin_52 = parse_tsplib_with_solution("berlin52.tsp", "berlin52.sol")
    dsj1000 = parse_tsplib_with_solution("dsj1000.tsp", "dsj1000.sol")
    # pla_85900 = parse_tsplib_with_solution("pla85900.tsp", "pla85900.sol")
    test_cases = [
        ("berlin_52", berlin_52),
        ("dsj1000", dsj1000),
        # ("pla_85900", None),
    ]

    idx = 0
    for name, (graph, expected_path, expected_cost) in test_cases:
        run_test_case(idx, graph, expected_path, expected_cost)
        idx += 1


if __name__ == '__main__':
    main_2()
