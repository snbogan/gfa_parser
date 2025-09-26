#!/usr/bin/env python3

### gfa_parser phased mode version 1.6 ###
### Compatible with hifiasm and verkko GFA files ###
### Written by Samuel N. Bogan [1], Owen W. Moosman [1], and Joanna L. Kelley [1] ###
### [1] University of California, Santa Cruz, Santa Cruz, USA ###

import argparse
import networkx as nx
import subprocess
import sys
import os


# Parse the GFA file into a directed graph and extract node/edge attributes
def parse_gfa(gfa_file):
    graph = nx.DiGraph()
    sequences, lengths, depths, orientations, edges = {}, {}, {}, {}, {}

    with open(gfa_file, 'r') as f:
        for line in f:
            # Segment line: defines a unitig and its properties
            if line.startswith('S'):
                parts = line.strip().split('\t')
                unitig_id = parts[1]
                sequence = parts[2]

                sequences[unitig_id] = sequence
                graph.add_node(unitig_id)

                # Extract length (LN) and read depth (rd), fallback if missing
                ln = next((int(p.split(':')[-1]) for p in parts if p.startswith('LN:i:')), len(sequence))
                rd = next((int(p.split(':')[-1]) for p in parts if p.startswith('rd:i:')), 0)

                lengths[unitig_id] = ln
                depths[unitig_id] = rd
                orientations[unitig_id] = '+'

            # Link line: defines an edge between two unitigs
            elif line.startswith('L'):
                parts = line.strip().split('\t')
                unitig1, strand1, unitig2, strand2, overlap = parts[1:6]

                # Extract overlap length (strip non-numeric chars)
                overlap_len = int(''.join(filter(str.isdigit, overlap)))

                graph.add_edge(unitig1, unitig2)
                edges[(unitig1, unitig2)] = (strand1, strand2, overlap_len)

    return graph, sequences, lengths, depths, orientations, edges



# Depth-first search to find all valid directed paths in the graph
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



### Build a full sequence from a path of unitigs
# Handles strand orientation and overlaps between connected unitigs
def build_sequence(path, sequences, edges):
    if not path:
        return ""

    current_seq = sequences[path[0]]
    for i in range(1, len(path)):
        prev, curr = path[i-1], path[i]
        strand1, strand2, overlap_len = edges.get((prev, curr), ('+', '+', 0))
        next_seq = sequences[curr]

        # Reverse complement if edge is on negative strand
        if strand1 == '-' or strand2 == '-':
            next_seq = next_seq[::-1].translate(str.maketrans("ATCGN", "TAGCN"))

        # Concatenate with overlap removed
        current_seq += next_seq[overlap_len:]

    return current_seq



### Construct haplotype-specific contigs:
# Filters unitigs below minimum read depth
# Finds valid paths from start to end unitigs
# Builds contig sequences and computes weighted read depth
def find_contigs(graph, sequences, lengths, depths, edges,
                 hap1_unitigs, hap2_unitigs,
                 hap1_start, hap1_end, hap2_start, hap2_end,
                 filter_rd):

    hap1_contigs, hap2_contigs = [], []

    # Filter unitigs by read depth threshold
    hap1_set = set(u for u in hap1_unitigs if depths.get(u, 0) >= filter_rd)
    hap2_set = set(u for u in hap2_unitigs if depths.get(u, 0) >= filter_rd)

    # Find haplotype 1 contigs
    hap1_paths = find_valid_paths(graph, hap1_start, hap1_end, hap1_set)
    for path in hap1_paths:
        if any(depths.get(u, 0) < filter_rd for u in path):
            continue
        seq = build_sequence(path, sequences, edges)
        total_len = sum(lengths[u] for u in path)
        weighted_rd = sum(lengths[u] * depths[u] for u in path) / total_len if total_len > 0 else 0
        hap1_contigs.append((path, seq, weighted_rd))

    # Find haplotype 2 contigs
    hap2_paths = find_valid_paths(graph, hap2_start, hap2_end, hap2_set)
    for path in hap2_paths:
        if any(depths.get(u, 0) < filter_rd for u in path):
            continue
        seq = build_sequence(path, sequences, edges)
        total_len = sum(lengths[u] for u in path)
        weighted_rd = sum(lengths[u] * depths[u] for u in path) / total_len if total_len > 0 else 0
        hap2_contigs.append((path, seq, weighted_rd))

    return hap1_contigs, hap2_contigs



# Write haplotype contigs to individual FASTA files
def write_individual_fasta(contigs, output_prefix, hap_label):
    for i, (path, seq, weighted_rd) in enumerate(contigs):
        contig_id = f"contig_{i+1}"
        output_file = f"{output_prefix}_hap{hap_label}_{contig_id}.fasta"
        with open(output_file, 'w') as f:
            f.write(f">{contig_id} weighted_rd={weighted_rd:.2f}\n")
            for j in range(0, len(seq), 60):  # wrap sequences at 60 bases per line
                f.write(seq[j:j+60] + '\n')



# Main program logic
def main(args):
    # If using shasta or minigraph, delegate to external scripts
    if args.package in ["shasta", "minigraph"]:
        script_name = "phased_mode_v1.6_shasta.py" if args.package == "shasta" else "phased_mode_v1.6_minigraph.py"
        script_path = os.path.join(os.path.dirname(__file__), script_name)
        cmd = [sys.executable, script_path] + sys.argv[1:]
        subprocess.run(cmd)
        return

    # Load haplotype-specific unitigs and start/end markers
    with open(args.hap1_unitigs) as f:
        hap1_unitigs = [line.strip() for line in f]
    with open(args.hap2_unitigs) as f:
        hap2_unitigs = [line.strip() for line in f]

    hap1_start = open(args.hap1_start).read().strip()
    hap1_end = open(args.hap1_end).read().strip()
    hap2_start = open(args.hap2_start).read().strip()
    hap2_end = open(args.hap2_end).read().strip()

    # Parse graph and reconstruct contigs
    graph, sequences, lengths, depths, orientations, edges = parse_gfa(args.gfa)
    hap1_contigs, hap2_contigs = find_contigs(
        graph, sequences, lengths, depths, edges,
        hap1_unitigs, hap2_unitigs,
        hap1_start, hap1_end,
        hap2_start, hap2_end,
        args.filter_rd
    )

    # Write haplotype-specific FASTA outputs
    write_individual_fasta(hap1_contigs, args.output_prefix, "1")
    write_individual_fasta(hap2_contigs, args.output_prefix, "2")



# Argument parsing
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract haplotype contigs from GFA paths.")

    parser.add_argument("--gfa", "-g", required=True, help="Input GFA file")
    parser.add_argument("--hap1-unitigs", "-h1u", required=True, help="File with unitig IDs for haplotype 1")
    parser.add_argument("--hap2-unitigs", "-h2u", required=True, help="File with unitig IDs for haplotype 2")
    parser.add_argument("--hap1-start", "-h1s", required=True, help="File with haplotype 1 start unitig ID")
    parser.add_argument("--hap1-end", "-h1e", required=True, help="File with haplotype 1 end unitig ID")
    parser.add_argument("--hap2-start", "-h2s", required=True, help="File with haplotype 2 start unitig ID")
    parser.add_argument("--hap2-end", "-h2e", required=True, help="File with haplotype 2 end unitig ID")
    parser.add_argument("--output-prefix", "-o", required=True, help="Prefix for output FASTA files")
    parser.add_argument("--filter_rd", "-r", type=int, default=0, help="Minimum read depth (rd) to keep unitigs (default: 0)")
    parser.add_argument("--package", "-p", choices=["hifiasm", "verkko", "shasta", "minigraph"], default="hifiasm",
                        help="Assembly package used (default: hifiasm)")

    args = parser.parse_args()
    main(args)



