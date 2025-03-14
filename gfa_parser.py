import sys
import networkx as nx
from collections import defaultdict

def parse_gfa(gfa_file):
    graph = nx.DiGraph()  # Use a directed graph
    sequences = {}
    with open(gfa_file, 'r') as f:
        for line in f:
            if line.startswith('S'):
                parts = line.strip().split('\t')
                unitig_id = parts[1]
                sequence = parts[2]
                sequences[unitig_id] = sequence
                graph.add_node(unitig_id)
            elif line.startswith('L'):
                parts = line.strip().split('\t')
                unitig1 = parts[1]
                unitig2 = parts[3]
                orientation1 = parts[2]
                orientation2 = parts[4]
                # Add directed edge based on the orientation
                graph.add_edge(unitig1, unitig2)
    return graph, sequences

def find_valid_paths(graph, start_unitig, end_unitig, valid_unitigs):
    """
    Find all valid directed paths from start_unitig to end_unitig
    that only include unitigs in valid_unitigs.
    """
    def dfs(current_node, path, visited):
        if current_node == end_unitig:
            valid_paths.append(path)
            return
        for neighbor in graph.successors(current_node):
            if neighbor in valid_unitigs and neighbor not in visited:
                dfs(neighbor, path + [neighbor], visited | {neighbor})

    valid_paths = []
    dfs(start_unitig, [start_unitig], {start_unitig})
    return valid_paths

def find_contigs(graph, sequences, hap1_unitigs, hap2_unitigs, hap1_start, hap1_end, hap2_start, hap2_end):
    hap1_contigs = []
    hap2_contigs = []
    
    # Create sets for quick lookup
    hap1_set = set(hap1_unitigs)
    hap2_set = set(hap2_unitigs)
    
    # Find valid paths for haplotype 1
    hap1_paths = find_valid_paths(graph, hap1_start, hap1_end, hap1_set)
    for path in hap1_paths:
        contig_seq = ''.join([sequences[node] for node in path])
        hap1_contigs.append((path, contig_seq))
    
    # Find valid paths for haplotype 2
    hap2_paths = find_valid_paths(graph, hap2_start, hap2_end, hap2_set)
    for path in hap2_paths:
        contig_seq = ''.join([sequences[node] for node in path])
        hap2_contigs.append((path, contig_seq))
    
    return hap1_contigs, hap2_contigs

def write_individual_fasta(contigs, output_prefix, hap_label):
    for i, (path, seq) in enumerate(contigs):
        contig_id = f"contig_{i+1}"
        output_file = f"{output_prefix}_hap{hap_label}_{contig_id}.fasta"
        with open(output_file, 'w') as f:
            f.write(f">{contig_id}\n")
            f.write(seq + "\n")

def main(gfa_file, hap1_unitigs_file, hap2_unitigs_file, hap1_start_file, hap1_end_file, hap2_start_file, hap2_end_file, output_prefix):
    with open(hap1_unitigs_file, 'r') as f:
        hap1_unitigs = [line.strip() for line in f]
    
    with open(hap2_unitigs_file, 'r') as f:
        hap2_unitigs = [line.strip() for line in f]

    with open(hap1_start_file, 'r') as f:
        hap1_start = f.readline().strip()
    
    with open(hap1_end_file, 'r') as f:
        hap1_end = f.readline().strip()
    
    with open(hap2_start_file, 'r') as f:
        hap2_start = f.readline().strip()
    
    with open(hap2_end_file, 'r') as f:
        hap2_end = f.readline().strip()

    graph, sequences = parse_gfa(gfa_file)
    hap1_contigs, hap2_contigs = find_contigs(graph, sequences, hap1_unitigs, hap2_unitigs, hap1_start, hap1_end, hap2_start, hap2_end)

    # Write individual FASTA files for each haplotype
    write_individual_fasta(hap1_contigs, output_prefix, "1")
    write_individual_fasta(hap2_contigs, output_prefix, "2")

if __name__ == "__main__":
    if len(sys.argv) != 9:
        print("Usage: python script.py <gfa_file> <hap1_unitigs_file> <hap2_unitigs_file> <hap1_start_file> <hap1_end_file> <hap2_start_file> <hap2_end_file> <output_prefix>")
        sys.exit(1)

    gfa_file = sys.argv[1]
    hap1_unitigs_file = sys.argv[2]
    hap2_unitigs_file = sys.argv[3]
    hap1_start_file = sys.argv[4]
    hap1_end_file = sys.argv[5]
    hap2_start_file = sys.argv[6]
    hap2_end_file = sys.argv[7]
    output_prefix = sys.argv[8]

    main(gfa_file, hap1_unitigs_file, hap2_unitigs_file, hap1_start_file, hap1_end_file, hap2_start_file, hap2_end_file, output_prefix)
