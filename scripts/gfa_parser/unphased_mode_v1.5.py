import argparse
import networkx as nx
import subprocess
import sys
import os

def parse_gfa(gfa_file):
    graph = nx.DiGraph()
    sequences = {}
    lengths = {}
    depths = {}
    orientations = {}  # Track original orientations
    edges = {}  # Store edge information including orientation
    
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
                orientations[unitig_id] = '+'  # Default orientation
            elif line.startswith('L'):
                parts = line.strip().split('\t')
                unitig1 = parts[1]
                strand1 = parts[2]
                unitig2 = parts[3]
                strand2 = parts[4]
                graph.add_edge(unitig1, unitig2)
                edges[(unitig1, unitig2)] = (strand1, strand2)
    
    return graph, sequences, lengths, depths, orientations, edges

# ... [keep all other functions the same: find_valid_paths, get_reverse_complement, build_sequence, find_contigs, write_individual_fasta] ...

def main(args):
    # Check package selection
    if args.package in ["shasta", "minigraph"]:
        script_to_run = "unphased_minigraph_v1.1.py" if args.package == "minigraph" else os.path.basename(__file__)
        cmd = [sys.executable, script_to_run] + sys.argv[1:]
        subprocess.run(cmd)
        return  # Exit after running the external script

    # Original script runs for hifiasm or verkko
    with open(args.hap1_unitigs, 'r') as f:
        hap1_unitigs = [line.strip() for line in f]

    with open(args.hap2_unitigs, 'r') as f:
        hap2_unitigs = [line.strip() for line in f]

    with open(args.hap1_start, 'r') as f:
        hap1_start = f.readline().strip()

    with open(args.hap1_end, 'r') as f:
        hap1_end = f.readline().strip()

    with open(args.hap2_start, 'r') as f:
        hap2_start = f.readline().strip()

    with open(args.hap2_end, 'r') as f:
        hap2_end = f.readline().strip()

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
    parser.add_argument("--filter_rd", "-f", type=int, default=0, help="Minimum rd to keep unitigs (default: 0)")
    parser.add_argument("--package", "-p", choices=["hifiasm", "verkko", "shasta", "minigraph"], default="hifiasm",
                        help="Which assembly package to run (default: hifiasm)")

    args = parser.parse_args()
    main(args)
