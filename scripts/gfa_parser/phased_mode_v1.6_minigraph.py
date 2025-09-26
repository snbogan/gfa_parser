#!/usr/bin/env python3

### gfa_parser phased mode version 1.6, compatible with minigraph gfa files ###
### Written by Samuel N. Bogan [1], Owen W. Moosman [1], and Joanna L. Kelley [1] ###
### [1] University of California, Santa Cruz, Santa Cruz, USA ###

import argparse
import sys
import os
import re
import networkx as nx
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# Parse GFA cigar strings
def parse_cigar_overlap(cigar):
    return sum(int(m) for m in re.findall(r'(\d+)M', cigar))

# Parse GFA networks
def parse_gfa_v3(gfa_file):
    sequences, lengths, depths, graph, edges = {}, {}, {}, nx.DiGraph(), {}
    with open(gfa_file,'r') as f:
        for line in f:
            if line.startswith('S'):
                parts = line.strip().split('\t')
                node = parts[1]; seq = parts[2]
                sequences[node]=seq; lengths[node]=len(seq)
                for tag in parts[3:]:
                    if tag.startswith('rd:i:'): depths[node]=int(tag.split(':')[-1])
            elif line.startswith('L'):
                parts=line.strip().split('\t')
                u1,s1,u2,s2,cigar=parts[1:6]
                ov_len=parse_cigar_overlap(cigar)
                graph.add_edge(u1,u2); edges[(u1,u2)] = (s1,s2,ov_len)
    return graph, sequences, lengths, depths, edges

# Build FASTAs of directed, acyclic paths accounting for overlaps and strandedness
def build_sequence(path,sequences,edges):
    if not path: return ""
    seq=sequences[path[0]]
    for i in range(1,len(path)):
        prev,curr=path[i-1],path[i]
        s1,s2,ov=edges.get((prev,curr),('+','+',0))
        next_seq=sequences[curr]
        if s2=='-': next_seq=str(Seq(next_seq).reverse_complement())
        seq+=next_seq[ov:]
    return seq

# Export FASTAs
def write_fastas(paths,sequences,depths,lengths,edges,filter_rd,outdir,hap_label):
    os.makedirs(outdir,exist_ok=True)
    count=0
    for path in paths:
        valid=[u for u in path if depths.get(u,0)>=filter_rd]
        if len(valid)<2: continue
        seq=build_sequence(valid,sequences,edges)
        total_len=sum(lengths[u] for u in valid)
        weighted_rd=sum(lengths[u]*depths[u] for u in valid)/total_len if total_len>0 else 0
        rec=SeqRecord(Seq(seq), id=f"contig_{hap_label}_{count}_rd{weighted_rd:.2f}", description="|".join(valid))
        with open(os.path.join(outdir,f"contig_{hap_label}_{count}.fasta"),'w') as handle: SeqIO.write(rec,handle,'fasta')
        count+=1

# Define arguments
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('-g','--gfa',required=True)
    parser.add_argument('-s','--start',required=True)
    parser.add_argument('-e','--end',required=True)
    parser.add_argument('-r','--filter_rd',type=int,default=0)
    parser.add_argument('-o','--outdir',required=True)
    args=parser.parse_args()
    start,end=open(args.start).read().strip(),open(args.end).read().strip()
    graph,sequences,lengths,depths,edges=parse_gfa_v3(args.gfa)
    paths=list(nx.all_simple_paths(graph,source=start,target=end))
    write_fastas(paths,sequences,depths,lengths,edges,args.filter_rd,args.outdir,"hap")

if __name__=='__main__':
    main()

