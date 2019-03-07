# TE_annotation
Transposable element annotation of Drosophila genomes

The first script that needs to be run is `annotate_TEs_polished.sh` .  

This script runs the top-level programs that find TEs in the given genome assembly: RepeatModeler, LTR finder, and dnaPipeTE. The way it is currently set up is specifically for the file names and directory structure of my current project - which will need to be modified in the extension of the pipeline. Currently the only input parameter it takes is the name of the species so it can find the file.  It is also set up to run on a Cornell CBSU-BSCB linux machine and the paths to the programs reflect this (all dependencies are already installed on the machine we used, and some programs are already in the machine's path). 

One thing that may need modification for larger genomes is the LTR finder module because it is very slow (~1 week) even for Drosophila-sized genomes. It is possible to break the genome into, for example, 20 Mb parts and send each part to a different node to be completed more efficiently. 

The second script that needs to be run after the first script is `TE_filter_cluster_wrapper.sh`.  
This script does the rest of the work in filtering the consensi, clustering them, and producing a final consensus library with minimum redundancy. This script uses helper scripts in the `src/` directory. 
Currently, this script takes in several input parameters: some of which can be removed to become automatic, and some input parameters that are built-in will need to be added to allow users more control. Currently it takes in the fileroot (everything but the .fasta) if the dnaPipeTE contigs, the RepeatModeler contigs, and the LTR_retriever sequences, and then the name of the species so it can find the appropriate folders and files.  

Current built-in parameters that should become user-specified parameters: max length of contig to be accepted for further processing: 10 kb, min length of contig to be accepted for further processing: 200 bp. The minimum cutoff is not applied to LTRs. These parameters were built based on Drosophila melanogaster and also all Repbase entries for other Drosophila species.  These parameters are used in the scripts: `src/filter_long_only.sh` and `src/filter_long_short.sh`.  

Another thing that will need to be modified for use in other genomes is the prefix used in the script `src/filter_LTR_clusters.sh`, specifically in the "for loop" on line 34. Here, in order to identify clusters that have LTR elements from LTR finder - I use the differntiating feature of LTR elements in that they are named based on their genomic coordinates, which begin with the prefix of the contig/scaffold names in the genome assembly. In my case, this prefix is "tig0". To make this applicable to other genomes, we will have to make a variable and determine the prefix(es) in the genome assembly empircally, or use a different method to identify LTR elements. 

The final output of the script is `consensi_clusters_final.filtered.fasta`.
