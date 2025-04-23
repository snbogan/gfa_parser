#!/usr/bin/env python3

import argparse
import sys
import os
import networkx as nx
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

def parse_gfa(gfa_file):
    sequences = {}
    lengths = {}
    depths = {}
    graph = nx.DiGraph()
    edges = {}

    with open(gfa_file, 'r') as f:
        for line in f:
            if line.startswith('S'):
                parts = line.strip().split('\t')
                node = parts[1]
                seq = parts[2]
                sequences[node] = seq
                lengths[node] = len(seq)
                for tag in parts[3:]:
                    if tag.startswith('rd:i:'):
                        try:
                            depths[node] = int(tag.split(':')[-1])
                        except ValueError:
                            depths[node] = 0
            elif line.startswith('L'):
                parts = line.strip().split('\t')
                from_node, from_strand, to_node, to_strand = parts[1:5]
                overlap = parts[5]  # e.g., 100M
                ov_len = int(overlap.rstrip('M'))
                graph.add_edge(from_node, to_node)
                edges[(from_node, to_node)] = (from_strand, to_strand, ov_len)

    return graph, sequences, lengths, depths, edges

def get_reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join([complement.get(base, 'N') for base in reversed(seq)])

def build_sequence(path, sequences, edges):
    if not path:
        return ""
    current_node = path[0]
    current_seq = sequences[current_node]

    for i in range(1, len(path)):
        prev_node = path[i-1]
        current_node = path[i]
        if (prev_node, current_node) in edges:
            strand1, strand2, overlap = edges[(prev_node, current_node)]
        elif (current_node, prev_node) in edges:
            strand2, strand1, overlap = edges[(current_node, prev_node)]
        else:
            strand1, strand2, overlap = '+', '+', 0  # assume no overlap

        next_seq = sequences[current_node]
        if strand2 == '-':
            next_seq = get_reverse_complement(next_seq)

        # Trim the overlap from the next sequence
        current_seq += next_seq[overlap:]

    return current_seq

def get_directed_paths(graph, start, end):
    if start not in graph or end not in graph:
        sys.exit(f"Error: Start or end node not found in graph ({start}, {end})")

    try:
        sp_len = nx.shortest_path_length(graph, source=start, target=end)
        paths = list(nx.all_simple_paths(graph, source=start, target=end, cutoff=sp_len))
    except nx.NetworkXNoPath:
        sys.exit("No directed path found between start and end.")

    if not paths:
        sys.exit("No directed shortest paths found.")

    unitig_set = set()
    for p in paths:
        unitig_set.update(p)
    unitig_set.discard(start)
    unitig_set.discard(end)

    return sorted(unitig_set), paths

def write_individual_fastas(paths, sequences, depths, lengths, edges, filter_rd, outdir):
    os.makedirs(outdir, exist_ok=True)
    count = 0

    for path in paths:
        valid_path = [u for u in path if depths.get(u, 0) >= filter_rd]
        if len(valid_path) < 2:
            continue
        seq = build_sequence(valid_path, sequences, edges)
        total_len = sum([lengths[u] for u in valid_path])
        weighted_rd = sum([lengths[u] * depths.get(u, 0) for u in valid_path]) / total_len if total_len > 0 else 0
        fasta_path = os.path.join(outdir, f'contig_{count}.fasta')
        record = SeqRecord(Seq(seq),
                           id=f"contig_{count}_rd{weighted_rd:.2f}",
                           description=f"path={'|'.join(valid_path)}")
        with open(fasta_path, 'w') as handle:
            SeqIO.write(record, handle, 'fasta')
        count += 1

    if count == 0:
        print("No contigs written (check read depth filtering).")
    else:
        print(f"Wrote {count} contig FASTA files to '{outdir}'.")

def main():
    parser = argparse.ArgumentParser(description="Extract directed paths between two flanking unitigs in a GFA and output individual FASTA files.")
    parser.add_argument('-g', '--gfa', required=True, help='Input GFA file')
    parser.add_argument('-s', '--start', required=True, help='File containing start unitig ID')
    parser.add_argument('-e', '--end', required=True, help='File containing end unitig ID')
    parser.add_argument('-r', '--filter_rd', type=int, default=0, help='Minimum read depth to include a unitig')
    parser.add_argument('-o', '--outdir', required=True, help='Directory to write output FASTA files')
    args = parser.parse_args()

    try:
        with open(args.start) as f:
            start = f.read().strip()
        with open(args.end) as f:
            end = f.read().strip()
    except Exception as e:
        sys.exit(f"Error reading start or end file: {e}")

    graph, sequences, lengths, depths, edges = parse_gfa(args.gfa)
    filtered_unitigs, paths = get_directed_paths(graph, start, end)
    write_individual_fastas(paths, sequences, depths, lengths, edges, args.filter_rd, args.outdir)

if __name__ == '__main__':
    main()
