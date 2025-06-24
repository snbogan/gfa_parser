# gfa_parser

In a genome assembly, a contig represents one assembly among multiple possible assemblies of shorter contiguous sequences called unitigs. Information about the possible networks of assemblies across unitigs are generally contained in graphical fragment assembly files (GFAs). 

*gfa_parser* is a tool to compute and extract all possible genome assemblies from GFA files using directed, acyclic networks.

## Authors 
Samuel N. Bogan | University of California, Santa Cruz

Owen W. Moosman | University of California, Santa Cruz

Joanna L. Kelley | University of California, Santa Cruz

## Method

Given a GFA file and set(s) of unitigs associated with a region or contig, *gfa_parser* identifies all possible contiguous paths through the flanking unitigs of each haplotype that are directed and acyclic using the [networkx package](https://networkx.org/). It also creates an index of unitig names and their associated sequences. The directed, acyclic paths ensure that (i) no assembly path loops through a unitig multiple times and (ii) that alternate unitigs within the same bubble are not included in the same path. The script then extracts the fasta sequence of each unitig and, for a given path, concatenates the fasta sequences in the order of the assembly path. This process is repeated for all possible assemblies of each haplotype. The script exports a directory of all possible contigs for both haplotypes in FASTA format. The FASTA header contains the length-weighted read depth of all unitigs included in the path and FASTA. A filtering parameter allows you to ignore unitigs and remove them from paths if they fall below a given read depth threshold.

The tool has two modes: unphased mode and phased mode. 

### Unphased mode

Computes and exports all directed acyclic paths through a gfa as fasta files, not separating based on haplotype.

### Phased mode

Computes and exports all directed acyclic paths through a gfa as fasta files for haplotypes 1 and 2 in a diploid genome.

## Unphased usage

    python unphased_mode.py \
    --gfa/-g my.gfa \ 			      # GFA file
    --start/-s start_unitig \ 		# ID of 5’ unitig
    --end/-e end_unitig \	      	# ID of 3’ unitig
    --filter_rd/-f 0 \			      # Skip unitigs below this read depth in paths
    --out/-o output_prefix 		    # Output
    
## Phased usage

    python phased_mode.py \
    --gfa/-g my.gfa \ 			            # GFA file
    --hap1_start/-h1s start_unitig \ 		# ID of hap 1 5’ unitig
    --hap1_end/-h1e end_unitig \	    	# ID of hap 1 3’ unitig
    --hap1_unitigs/-h1u \ 				      # List of hap 1 unitig IDs to consider
    --hap2_start/-h2s start_unitig \ 		# ID of hap 1 5’ unitig
    --hap2_end/-h2e end_unitig \		    # ID of hap 2 3’ unitig
    --hap2_unitigs/-h2u \ 				      # List of hap 2 unitig IDs to consider
    --filter_rd/-f 0 \			            # Skip unitigs below this read depth in paths
    --out/-o output_prefix 		          # Output
    
## Citation

A manuscript reporting *gfa_parser* is in review. For now, please cite this Github page and the listed authors if you use *gfa_parser*.
    
