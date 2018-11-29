#! /bin/bash

# filtering and clustering TE annotation wrapper script
# use the DPTE filtered (for genomic abundance) output, and Rmod sequences (not the annotated ones)

# TE_filter_cluster_wrapper.sh <DPTE_seqs> <RMod_seqs> <LTR_retriever> <species>

echo "dnaPipeTE sequences:"
echo $1
echo "Repeat Modeler sequences:"
echo $2
echo "LTR retriever sequences:"
echo $3
echo "species:"
echo $4

# first, filter sequences that are too long or too short

mkdir extra_files

# for DPTE and RMod - remove both long and short
# for LTR retriever, remove only too long

sh src/filter_long_short.sh $1 200
sh src/filter_long_short.sh $2 200
sh src/filter_long_only.sh $3 

mv *long extra_files
mv *short extra_files
mv *badlength* extra_files

# next, remove the sequences that have Ns from the Rmod output

python src/removeNfromfas.py $2.goodlength.fasta
mv N_removed.fasta $2.noNs.fasta

### TAKE CARE OF THE LTRs by LTR_finder/LTR_retriever  ###

# run mothur and get OTU representative
sh src/mothur_LTRs.sh $3.goodlength
# output file is $3.LTR.OTUreps.fasta
mv *lengths *candidates* *singleton* extra_files

### Now cluster the LTR OTUs with the filtered DPTE and Rmod output, then remove the ones that cluster with these OTUs, and then re-cluster #####

sh src/filter_LTR_clusters.sh $1.goodlength $2.noNs $3.goodlength.LTR.OTUreps  # this script will have to be modified for Dmel vs. Dvirilis. 

#output is consensi_clusters_final.fasta

# move the extra output files to extra folder
mv seqs* *clstr extra_files

# done.





