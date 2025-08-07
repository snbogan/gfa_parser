#!/usr/bin/env python3

import networkx as nx
import argparse
from tqdm import tqdm
from decimal import Decimal, getcontext
import sys
import math

def parse_gfa(gfa_file):
    graph = nx.DiGraph()
    lengths = {}

    with open(gfa_file, 'r') as f:
        for line in f:
            if line.startswith('S'):
                parts = line.strip().split('\t')
                node = parts[1]
                seq = parts[2]
                graph.add_node(node)
                lengths[node] = len(seq)
            elif line.startswith('L'):
                parts = line.strip().split('\t')
                from_node = parts[1]
                to_node = parts[3]
                graph.add_edge(from_node, to_node)

    return graph, lengths

def count_paths_from_node(graph, start_node, memo):
    if start_node in memo:
        return memo[start_node]

    stack = [(start_node, iter(graph.successors(start_node)))]
    path_counts = {start_node: 0}

    while stack:
        current, children = stack[-1]
        try:
            child = next(children)
            if child in memo:
                path_counts[current] += 1 + memo[child]
            elif child not in path_counts:
                path_counts[child] = 0
                stack.append((child, iter(graph.successors(child))))
        except StopIteration:
            stack.pop()
            for neighbor in graph.successors(current):
                path_counts[current] += 1 + path_counts.get(neighbor, 0)

    memo[start_node] = path_counts[start_node]
    return path_counts[start_node]

def count_all_dag_paths(graph):
    total_count = 0
    memo = {}

    for node in tqdm(graph.nodes(), desc="Traversing unitigs", unit="unitig"):
        total_count += count_paths_from_node(graph, node, memo)

    return total_count

# Expected log10(DAGs) from number of unitigs, using fitted or assumed model
def expected_log10_dags(num_unitigs):
    # Example model coefficients; update these with empirical fitting
    a = 0.015
    b = -0.1
    c = 2.0
    if num_unitigs <= 1:
        return 0
    return a * num_unitigs * math.log(num_unitigs) + b * num_unitigs + c

def main():
    parser = argparse.ArgumentParser(description="Efficient count of DAG paths in a GFA using overlaps.")
    parser.add_argument('-g', '--gfa', required=True, help='Input GFA file')
    args = parser.parse_args()

    graph, lengths = parse_gfa(args.gfa)
    total_length = sum(lengths.values())
    num_unitigs = len(lengths)

    print("Counting DAG paths only through reachable overlaps...")
    total_paths = count_all_dag_paths(graph)

    getcontext().prec = 50  # Prevent float overflow with very large integers
    total_paths_decimal = Decimal(total_paths)
    num_unitigs_decimal = Decimal(num_unitigs)

    print(f"\nTotal unitigs: {num_unitigs}")
    print(f"Total unitig length: {total_length}")
    print("Total directed acyclic paths:", f"{total_paths_decimal:.3e}")

    normalized = total_paths_decimal / Decimal(total_length)
    print(f"Normalized path count (paths per base): {normalized:.3e}")

    avg_paths_per_unitig = total_paths_decimal / num_unitigs_decimal
    print(f"Average number of DAGs per unitig: {avg_paths_per_unitig:.3e}")

    # Observed log10(DAGs)
    log10_dags = math.log10(total_paths)
    print(f"log10(total DAGs): {log10_dags:.3f}")

    # Expected log10(DAGs) from model
    expected_log_dags = expected_log10_dags(num_unitigs)
    print(f"Expected log10(DAGs) from model: {expected_log_dags:.3f}")

    # Residual = observed - expected
    residual = log10_dags - expected_log_dags
    print(f"Normalized DAG residual (obs - exp): {residual:.3f}")

if __name__ == '__main__':
    sys.setrecursionlimit(100000)
    main()
