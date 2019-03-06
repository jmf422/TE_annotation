#! /bin/bash

# filtering and clustering TE annotation wrapper script
# use the DPTE filtered (for genomic abundance) output, and Rmod sequences (not the annotated ones)

# TE_filter_cluster_wrapper.sh <DPTE_seqs> <RMod_seqs> <LTR_retriever> <species>
# consensi.fasta D_virilis.DPTE.contigs D_virilis.LTRlib

echo "dnaPipeTE sequences:"
echo $1
echo "Repeat Modeler sequences:"
echo $2
echo "LTR retriever sequences:"
echo $3
echo "species:"
echo $4


# first, filter sequences that are too long or too short
src="/workdir/jmf422/Polished_genomes/src"

cd $4


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

# next, get the sequences with Ns

python $src/removeNfromfas.py $2.goodlength.fasta
mv N_removed.fasta $2.noNs.fasta
cat $2.goodlength.fasta | grep '^>' | cut -f 2 -d '>' > $2.goodlength
cat $2.noNs.fasta | grep '^>' | cut -f 2 -d '>' > $2.noNs

sort $2.goodlength > $2.goodlength.sorted
sort $2.noNs > $2.noNs.sorted
comm -23 $2.goodlength.sorted $2.noNs.sorted > $2.Ns
echo "got Rmod sequences with Ns"
xargs samtools faidx $2.goodlength.fasta < $2.Ns > $2.Ns.fasta

echo "got fasta file of sequences with Ns"

# then re-assemble the sequences with Ns
sh $src/assemble_seqs_Ns_pipeline.sh $2.Ns $4.assembly $4

echo "done assembly of sequences with Ns"

#output is Nseq.contigs.passed.fasta, already filtered

# combine the seqs with Ns to the seqs without Ns
cat $2.noNs.fasta Nseq.contigs.passed.fasta > $2.final.fasta

cp $2.final.fasta ../final_files


### TAKE CARE OF THE LTRs by LTR_finder/LTR_retriever  ###
echo "now taking care of LTR library"
# run mothur and get OTU representative
sh $src/mothur_LTRs.sh $3.goodlength
# output file is $3.LTR.OTUreps.fasta

cp $3.LTR.OTUreps.fasta ../final_files

### Now cluster the LTR OTUs with the filtered DPTE and Rmod output, then remove the ones that cluster with these OTUs, and then re-cluster #####

sh $src/filter_LTR_clusters.sh $1.goodlength $2.final $3.goodlength.LTR.OTUreps  # this script will have to be modified for Dmel vs. Dvirilis. 

#output is consensi_clusters_final.fasta

# now filter the consensi based on if they mask the genome at least 1 time.
echo "now filter the consensi based on if they mask the genome at least 1 time"

sh $src/eval_consensi_Rmasker.sh consensi_clusters_final $4.assembly $4

# final output is consensi_clusters_final.filtered.fasta

mv consensi_clusters_final.filtered.fasta ../final_files




