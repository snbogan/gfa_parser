import argparse
import networkx as nx

def parse_gfa(gfa_file):
    graph = nx.DiGraph()
    sequences = {}
    lengths = {}
    depths = {}
    orientations = {}  # Track original orientations
    edges = {}  # Store edge information including orientation and overlap

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
                overlap = parts[5]  # Format: e.g., 20M
                overlap_len = int(''.join(filter(str.isdigit, overlap)))  # extract numeric part
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

def get_reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join([complement.get(base, 'N') for base in reversed(seq)])

def build_sequence(path, sequences, edges):
    if not path:
        return ""
    
    # Initialize with the first node
    first_node = path[0]
    current_seq = sequences[first_node]
    orientation = '+'

    for i in range(1, len(path)):
        prev_node = path[i-1]
        current_node = path[i]

        strand1, strand2, overlap_len = edges.get((prev_node, current_node), ('+', '+', 0))

        # Handle reverse complements
        prev_seq = current_seq
        next_seq = sequences[current_node]
        if strand1 == '-' and strand2 == '-':
            next_seq = get_reverse_complement(next_seq)
        elif strand1 == '+' and strand2 == '-':
            next_seq = get_reverse_complement(next_seq)
        elif strand1 == '-' and strand2 == '+':
            next_seq = get_reverse_complement(next_seq)
        # otherwise both '+' so use as is

        # Trim overlap from the start of the current unitig
        current_seq += next_seq[overlap_len:]

    return current_seq

def find_contigs(graph, sequences, lengths, depths, edges,
                 hap1_unitigs, hap2_unitigs,
                 hap1_start, hap1_end, hap2_start, hap2_end,
                 filter_rd):
    hap1_contigs = []
    hap2_contigs = []

    hap1_set = set(u for u in hap1_unitigs if depths.get(u, 0) >= filter_rd)
    hap2_set = set(u for u in hap2_unitigs if depths.get(u, 0) >= filter_rd)

    hap1_paths = find_valid_paths(graph, hap1_start, hap1_end, hap1_set)
    for path in hap1_paths:
        if any(depths.get(u, 0) < filter_rd for u in path):
            continue
        contig_seq = build_sequence(path, sequences, edges)
        total_len = sum([lengths[node] for node in path])
        weighted_rd = sum([lengths[node] * depths[node] for node in path]) / total_len if total_len > 0 else 0
        hap1_contigs.append((path, contig_seq, weighted_rd))

    hap2_paths = find_valid_paths(graph, hap2_start, hap2_end, hap2_set)
    for path in hap2_paths:
        if any(depths.get(u, 0) < filter_rd for u in path):
            continue
        contig_seq = build_sequence(path, sequences, edges)
        total_len = sum([lengths[node] for node in path])
        weighted_rd = sum([lengths[node] * depths[node] for node in path]) / total_len if total_len > 0 else 0
        hap2_contigs.append((path, contig_seq, weighted_rd))

    return hap1_contigs, hap2_contigs

def write_individual_fasta(contigs, output_prefix, hap_label):
    for i, (path, seq, weighted_rd) in enumerate(contigs):
        contig_id = f"contig_{i+1}"
        output_file = f"{output_prefix}_hap{hap_label}_{contig_id}.fasta"
        with open(output_file, 'w') as f:
            f.write(f">{contig_id} weighted_rd={weighted_rd:.2f}\n")
            for j in range(0, len(seq), 60):
                f.write(seq[j:j+60] + '\n')

def main(args):
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

    args = parser.parse_args()
    main(args)
