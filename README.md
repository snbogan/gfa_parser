# gfa_parser

In a genome, a contig represents one assembly among multiple possible assemblies of shorter contiguous sequences called unitigs. Information about the possible networks of assemblies across unitigs are generally contained in graphical fragment assembly files (GFAs). 

*gfa_parser* is a tool to compute and extract all possible genome assemblies from GFA files using directed, acyclic networks.

## Authors: 
Samuel N. Bogan | University of California, Santa Cruz

Owen W. Moosman | University of California, Santa Cruz

Joanna L. Kelley | University of California, Santa Cruz

## Method

Given a GFA file and set(s) of unitigs associated with a region or contig, *gfa_parser* identifies all possible contiguous paths through the flanking unitigs of each haplotype that are directed and acyclic using the [networkx package](https://networkx.org/). It also creates an index of unitig names and their associated sequences. The directed, acyclic paths ensure that (i) no assembly path loops through a unitig multiple times and (ii) that alternate unitigs within the same bubble are not included in the same path. The script then extracts the fasta sequence of each unitig and, for a given path, concatenates the fasta sequences in the order of the assembly path. This process is repeated for all possible assemblies of each haplotype. The script exports a directory of all possible contigs for both haplotypes in FASTA format.

## Usage

The script gfa_parser.py reads in four sets of arguments: 

1. A GFA file contain unitig-associated FASTA sequences 

2. Two .txt lists of unitigs (utg) phased to each haplotype in a contig of interest

3.. ‘start’ and ‘end’ unitigs flanking the 5’ and 3’ ends of each haplotype

4. A prefix for output FASTA files. 

It can be run as follows:

    python gfa_parser.py [gfa file] \
    [hap1 unitigs] [hap2 utg's] \
    [hap1 5' utg's] [hap1 3' utg's] \
    [hap2 5' utg's] [hap2 3' utg's] \
    [output prefix]
    
## Citation

A manuscript reporting *gfa_parser* is in review. For now, please cite this Github page and the listed authors if you use *gfa_parser*.
    
