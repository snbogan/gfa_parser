#!/usr/bin/env python3
import argparse
import networkx as nx
import subprocess
import sys
import os

# ----------------------------
# Original Script 4 logic (hifiasm/verkko)
# ----------------------------
def parse_gfa(gfa_file):
    graph = nx.DiGraph()
    sequences, lengths, depths, orientations, edges = {}, {}, {}, {}, {}
    with open(gfa_file, 'r') as f:
        for line in f:
            if line.startswith('S'):
                parts = line.strip().split('\t')
                unitig_id = parts[1]
                sequence = parts[2]
                sequences[unitig_id] = sequence
                graph.add_node(unitig_id)
                ln = next((int(p.split(':')[-1]) for p in parts if p.startswith('LN:i:')), len(sequence))
                rd = next((int(p.split(':')[-1]) for p in parts if p.startswith('rd:i:')), 0)
                lengths[unitig_id] = ln
                depths[unitig_id] = rd
                orientations[unitig_id] = '+'
            elif line.startswith('L'):
                parts = line.strip().split('\t')
                unitig1, strand1, unitig2, strand2, overlap = parts[1:6]
                overlap_len = int(''.join(filter(str.isdigit, overlap)))
                graph.add_edge(unitig1, unitig2)
                edges[(unitig1, unitig2)] = (strand1, strand2, overlap_len)
    return graph, sequences, lengths, depths, orientations, edges

def find_valid_paths(graph, start_unitig, end_unitig, valid_unitigs):
    def dfs(current_node, path, visited):
        if current_node == end_unitig:
            valid_paths.append(list(path))
            return
        for neighbor in graph.successors(current_node):
            if neighbor in valid_unitigs and neighbor not in visited:
                visited.add(neighbor)
                path.append(neighbor)
                dfs(neighbor, path, visited)
                path.pop()
                visited.remove(neighbor)

    valid_paths = []
    dfs(start_unitig, [start_unitig], {start_unitig})
    return valid_paths

def build_sequence(path, sequences, edges):
    if not path:
        return ""
    current_seq = sequences[path[0]]
    for i in range(1, len(path)):
        prev, curr = path[i-1], path[i]
        strand1, strand2, overlap_len = edges.get((prev, curr), ('+', '+', 0))
        next_seq = sequences[curr]
        if strand1 == '-' or strand2 == '-':
            next_seq = next_seq[::-1].translate(str.maketrans("ATCGN", "TAGCN"))
        current_seq += next_seq[overlap_len:]
    return current_seq

def find_contigs(graph, sequences, lengths, depths, edges,
                 hap1_unitigs, hap2_unitigs,
                 hap1_start, hap1_end, hap2_start, hap2_end,
                 filter_rd):
    hap1_contigs, hap2_contigs = [], []
    hap1_set = set(u for u in hap1_unitigs if depths.get(u, 0) >= filter_rd)
    hap2_set = set(u for u in hap2_unitigs if depths.get(u, 0) >= filter_rd)

    hap1_paths = find_valid_paths(graph, hap1_start, hap1_end, hap1_set)
    for path in hap1_paths:
        if any(depths.get(u, 0) < filter_rd for u in path):
            continue
        seq = build_sequence(path, sequences, edges)
        total_len = sum(lengths[u] for u in path)
        weighted_rd = sum(lengths[u] * depths[u] for u in path) / total_len if total_len > 0 else 0
        hap1_contigs.append((path, seq, weighted_rd))

    hap2_paths = find_valid_paths(graph, hap2_start, hap2_end, hap2_set)
    for path in hap2_paths:
        if any(depths.get(u, 0) < filter_rd for u in path):
            continue
        seq = build_sequence(path, sequences, edges)
        total_len = sum(lengths[u] for u in path)
        weighted_rd = sum(lengths[u] * depths[u] for u in path) / total_len if total_len > 0 else 0
        hap2_contigs.append((path, seq, weighted_rd))

    return hap1_contigs, hap2_contigs

def write_individual_fasta(contigs, output_prefix, hap_label):
    for i, (path, seq, weighted_rd) in enumerate(contigs):
        contig_id = f"contig_{i+1}"
        output_file = f"{output_prefix}_hap{hap_label}_{contig_id}.fasta"
        with open(output_file, 'w') as f:
            f.write(f">{contig_id} weighted_rd={weighted_rd:.2f}\n")
            for j in range(0, len(seq), 60):
                f.write(seq[j:j+60] + '\n')

# ----------------------------
# Main
# ----------------------------
def main(args):
    # If user selects shasta or minigraph, call the external script
    if args.package in ["shasta", "minigraph"]:
        script_name = "phased_mode_v1.6_shasta.py" if args.package == "shasta" else "phased_mode_v1.6_minigraph.py"
        script_path = os.path.join(os.path.dirname(__file__), script_name)
        cmd = [sys.executable, script_path] + sys.argv[1:]
        subprocess.run(cmd)
        return

    # Original hifiasm/verkko functionality
    with open(args.hap1_unitigs) as f:
        hap1_unitigs = [line.strip() for line in f]
    with open(args.hap2_unitigs) as f:
        hap2_unitigs = [line.strip() for line in f]

    hap1_start = open(args.hap1_start).read().strip()
    hap1_end = open(args.hap1_end).read().strip()
    hap2_start = open(args.hap2_start).read().strip()
    hap2_end = open(args.hap2_end).read().strip()

    graph, sequences, lengths, depths, orientations, edges = parse_gfa(args.gfa)
    hap1_contigs, hap2_contigs = find_contigs(
        graph, sequences, lengths, depths, edges,
        hap1_unitigs, hap2_unitigs,
        hap1_start, hap1_end,
        hap2_start, hap2_end,
        args.filter_rd
    )

    write_individual_fasta(hap1_contigs, args.output_prefix, "1")
    write_individual_fasta(hap2_contigs, args.output_prefix, "2")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract contigs from GFA based on haplotype paths.")
    parser.add_argument("--gfa", "-g", required=True, help="Input GFA file")
    parser.add_argument("--hap1-unitigs", "-h1u", required=True, help="File with list of unitigs for haplotype 1")
    parser.add_argument("--hap2-unitigs", "-h2u", required=True, help="File with list of unitigs for haplotype 2")
    parser.add_argument("--hap1-start", "-h1s", required=True, help="File containing haplotype 1 start unitig")
    parser.add_argument("--hap1-end", "-h1e", required=True, help="File containing haplotype 1 end unitig")
    parser.add_argument("--hap2-start", "-h2s", required=True, help="File containing haplotype 2 start unitig")
    parser.add_argument("--hap2-end", "-h2e", required=True, help="File containing haplotype 2 end unitig")
    parser.add_argument("--output-prefix", "-o", required=True, help="Prefix for output FASTA files")
    parser.add_argument("--filter_rd", "-r", type=int, default=0, help="Minimum rd to keep unitigs (default: 0)")
    parser.add_argument("--package", "-p", choices=["hifiasm", "verkko", "shasta", "minigraph"], default="hifiasm",
                        help="Which assembly package to run (default: hifiasm)")

    args = parser.parse_args()
    main(args)


