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
    depths = {}       # read depth (hifiasm rd:i: or fallback)
    sr_counts = {}    # Shasta SR:i: counts for each unitig
    graph = nx.DiGraph()
    edges = {}

    with open(gfa_file, 'r') as f:
        for line in f:
            if line.startswith('S'):  # Segment line
                parts = line.strip().split('\t')
                node = parts[1]
                seq = parts[2]
                sequences[node] = seq
                lengths[node] = len(seq)

                depth_val = None
                sr_val = None

                for tag in parts[3:]:
                    if tag.startswith('rd:i:'):  # hifiasm read depth
                        try:
                            depth_val = int(tag.split(':')[-1])
                        except ValueError:
                            depth_val = 0
                    elif tag.startswith('SR:i:'):  # Shasta read count
                        sr_val = int(tag.split(':')[-1])

                # Fallbacks if rd:i: missing
                if depth_val is None:
                    if sr_val is not None:
                        depth_val = sr_val  # use SR count as coverage proxy
                    else:
                        depth_val = 0

                depths[node] = depth_val
                sr_counts[node] = sr_val if sr_val is not None else 0

            elif line.startswith('L'):  # Link line
                parts = line.strip().split('\t')
                from_node, from_strand, to_node, to_strand = parts[1:5]
                overlap = parts[5]  # e.g., 100M or 0M
                try:
                    ov_len = int(overlap.rstrip('M'))
                except ValueError:
                    ov_len = 0
                graph.add_edge(from_node, to_node)
                edges[(from_node, to_node)] = (from_strand, to_strand, ov_len)

    return graph, sequences, lengths, depths, sr_counts, edges

def get_reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join([complement.get(base, 'N') for base in reversed(seq)])

def build_sequence(path, sequences, edges):
    if not path:
        return ""
    current_node = path[0]
    current_seq = sequences[current_node]

    for i in range(1, len(path)):
        prev_node = path[i - 1]
        current_node = path[i]
        if (prev_node, current_node) in edges:
            strand1, strand2, overlap = edges[(prev_node, current_node)]
        elif (current_node, prev_node) in edges:
            strand2, strand1, overlap = edges[(current_node, prev_node)]
        else:
            strand1, strand2, overlap = '+', '+', 0

        next_seq = sequences[current_node]
        if strand2 == '-':
            next_seq = get_reverse_complement(next_seq)

        current_seq += next_seq[overlap:] if overlap < len(next_seq) else ""

    return current_seq

def get_all_paths(graph, start, end, max_paths=None):
    if start not in graph or end not in graph:
        sys.exit(f"Error: Start or end node not found in graph ({start}, {end})")

    try:
        paths = nx.all_simple_paths(graph, source=start, target=end)
        if max_paths:
            paths = (p for i, p in enumerate(paths) if i < max_paths)
        else:
            paths = list(paths)
    except nx.NetworkXNoPath:
        sys.exit("No directed path found between start and end.")

    return paths

def write_individual_fastas(paths, sequences, depths, sr_counts, lengths, edges, filter_rd, outdir):
    os.makedirs(outdir, exist_ok=True)
    count = 0

    for path in paths:
        valid_path = [u for u in path if depths.get(u, 0) >= filter_rd]
        if len(valid_path) < 2:
            continue

        seq = build_sequence(valid_path, sequences, edges)
        total_len = sum([lengths[u] for u in valid_path])

        weighted_rd = sum([lengths[u] * depths.get(u, 0) for u in valid_path]) / total_len if total_len > 0 else 0
        mean_sr = sum([sr_counts.get(u, 0) for u in valid_path]) / len(valid_path) if valid_path else 0

        fasta_path = os.path.join(outdir, f'contig_{count}.fasta')
        record = SeqRecord(
            Seq(seq),
            id=f"contig_{count}_rd{weighted_rd:.2f}_meanSR{mean_sr:.2f}",
            description=f"path={'|'.join(valid_path)}"
        )
        with open(fasta_path, 'w') as handle:
            SeqIO.write(record, handle, 'fasta')
        count += 1

    if count == 0:
        print("No contigs written (check read depth filtering).")
    else:
        print(f"Wrote {count} contig FASTA files to '{outdir}'.")

def main():
    parser = argparse.ArgumentParser(description="Extract all directed paths between two flanking unitigs in a GFA and output individual FASTA files.")
    parser.add_argument('-g', '--gfa', required=True, help='Input GFA file')
    parser.add_argument('-s', '--start', required=True, help='File containing start unitig ID')
    parser.add_argument('-e', '--end', required=True, help='File containing end unitig ID')
    parser.add_argument('-r', '--filter_rd', type=int, default=0, help='Minimum read depth to include a unitig')
    parser.add_argument('-o', '--outdir', required=True, help='Directory to write output FASTA files')
    parser.add_argument('--max_paths', type=int, default=None, help='Optional: limit number of paths to enumerate (avoid huge memory usage)')
    args = parser.parse_args()

    try:
        with open(args.start) as f:
            start = f.read().strip()
        with open(args.end) as f:
            end = f.read().strip()
    except Exception as e:
        sys.exit(f"Error reading start or end file: {e}")

    graph, sequences, lengths, depths, sr_counts, edges = parse_gfa(args.gfa)
    paths = get_all_paths(graph, start, end, max_paths=args.max_paths)
    write_individual_fastas(paths, sequences, depths, sr_counts, lengths, edges, args.filter_rd, args.outdir)

if __name__ == '__main__':
    main()

