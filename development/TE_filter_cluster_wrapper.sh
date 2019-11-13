#! /bin/bash

# filtering and clustering TE annotation wrapper script
# use the DPTE filtered (for genomic abundance) output, and Rmod sequences (not the annotated ones)

# do from /workdir/jmf422/Final_annotations/Dmel_6-16-2019

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
src="/workdir/jmf422/Final_annotations/Dmel_6-16-2019/src"

#cd $4


mkdir clustering_annotation
mkdir final_files
cd clustering_annotation

cp ../Annotation_files/$1.fasta .
cp ../Annotation_files/$2.fasta .
cp ../Annotation_files/$3.fasta .


# for DPTE and RMod - remove both long and short
# for LTR retriever, remove only too long

sh $src/filter_long_short.sh $1 200
sh $src/filter_long_short.sh $2 200
sh $src/filter_long_only.sh $3 


echo "done first length filtering" 


# $2.goodlength.fasta - new LTR retreiver



### TAKE CARE OF THE LTRs by LTR_finder/LTR_retriever  ###
echo "now taking care of LTR library"
# run mothur and get OTU representative
sh $src/mothur_LTRs.sh $3.goodlength
#output: clustered_refined_LTRlib.fasta

cp clustered_refined_LTRlib.fasta ../final_files

### Now cluster the LTR OTUs with the filtered DPTE and Rmod output, then remove the ones that cluster with these OTUs, and then re-cluster #####

sh $src/filter_LTR_clusters.sh $1.goodlength $2.goodlength clustered_refined_LTRlib

#output is consensi_clusters_final.fasta



# now filter the consensi based on if they mask the genome at least 1 time.
echo "now filter the consensi based on if they mask the genome at least 1 times"

sh $src/eval_consensi_Rmasker.sh consensi_clusters_final $4.assembly $4

# final output is consensi_clusters_final.filtered.fasta

mv consensi_clusters_final.filtered.fasta ../final_files