#!/usr/bin/env python3
import argparse
import sys
import os
import networkx as nx
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def parse_gfa_v2(gfa_file):
    sequences, lengths, depths, sr_counts, graph, edges = {}, {}, {}, {}, nx.DiGraph(), {}
    with open(gfa_file, 'r') as f:
        for line in f:
            if line.startswith('S'):
                parts = line.strip().split('\t')
                node = parts[1]
                seq = parts[2]
                sequences[node] = seq
                lengths[node] = len(seq)
                depth_val, sr_val = None, None
                for tag in parts[3:]:
                    if tag.startswith('rd:i:'):
                        depth_val = int(tag.split(':')[-1])
                    elif tag.startswith('SR:i:'):
                        sr_val = int(tag.split(':')[-1])
                depths[node] = depth_val if depth_val is not None else (sr_val if sr_val is not None else 0)
                sr_counts[node] = sr_val if sr_val is not None else 0
            elif line.startswith('L'):
                parts = line.strip().split('\t')
                u1, s1, u2, s2, overlap = parts[1:6]
                ov_len = int(overlap.rstrip('M')) if overlap.endswith('M') else 0
                graph.add_edge(u1, u2)
                edges[(u1, u2)] = (s1, s2, ov_len)
    return graph, sequences, lengths, depths, sr_counts, edges

def build_sequence(path, sequences, edges):
    if not path: return ""
    seq = sequences[path[0]]
    for i in range(1, len(path)):
        prev, curr = path[i-1], path[i]
        s1, s2, ov = edges.get((prev, curr), ('+', '+', 0))
        next_seq = sequences[curr]
        if s2 == '-': next_seq = str(Seq(next_seq).reverse_complement())
        seq += next_seq[ov:]
    return seq

def write_fastas(paths, sequences, depths, sr_counts, lengths, edges, filter_rd, outdir, hap_label):
    os.makedirs(outdir, exist_ok=True)
    count = 0
    for path in paths:
        valid_path = [u for u in path if depths.get(u,0) >= filter_rd]
        if len(valid_path) < 2: continue
        seq = build_sequence(valid_path, sequences, edges)
        total_len = sum(lengths[u] for u in valid_path)
        weighted_rd = sum(lengths[u]*depths.get(u,0) for u in valid_path)/total_len if total_len>0 else 0
        mean_sr = sum(sr_counts.get(u,0) for u in valid_path)/len(valid_path)
        path_id = f"{hap_label}_{count}"
        rec = SeqRecord(Seq(seq), id=f"contig_{path_id}_rd{weighted_rd:.2f}_meanSR{mean_sr:.2f}",
                        description="|".join(valid_path))
        with open(os.path.join(outdir,f"contig_{path_id}.fasta"),'w') as handle: SeqIO.write(rec, handle,'fasta')
        count += 1

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g','--gfa', required=True)
    parser.add_argument('-s','--start', required=True)
    parser.add_argument('-e','--end', required=True)
    parser.add_argument('-r','--filter_rd', type=int, default=0)
    parser.add_argument('-o','--outdir', required=True)
    args = parser.parse_args()
    start, end = open(args.start).read().strip(), open(args.end).read().strip()
    graph, sequences, lengths, depths, sr_counts, edges = parse_gfa_v2(args.gfa)
    import networkx as nx
    paths = list(nx.all_simple_paths(graph, source=start, target=end))
    write_fastas(paths, sequences, depths, sr_counts, lengths, edges, args.filter_rd, args.outdir, "hap")

if __name__=='__main__':
    main()

