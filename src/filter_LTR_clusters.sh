#! /bin/bash
# Remove the redundant LTR sequences between LTR retriever and the other programs, and then re-cluster.

# sh filter_LTR_clusters.sh <DPTE sequences> <Rmod sequences> <LTR OTU reps>

# we need to modify the names of DPTE
cat $1.fasta | sed 's/comp_TRINITY_//g'  > temp
mv temp $1.fasta


# first combine the files
cat $1.fasta $2.fasta $3.fasta > all_consensi.fasta

# first run CDhit
/programs/cd-hit-v4.6.1-2012-08-27/cd-hit-est -aS 0.8 -c 0.8 -g 1 -G 0 -A 80 -M 10000 -i all_consensi.fasta -o all_consensi.clusters1.fasta -T 20

# then break up the file to process it.

csplit --digits=4 --quiet --prefix=cluster  all_consensi.clusters1.fasta.clstr "/^>/" "{*}"

# 718000 is common to all LTRs in this case

# in the case of the merged assemblies, it will be either tig or 000

# grep return T or F
#https://unix.stackexchange.com/questions/48535/can-grep-return-true-false-or-are-there-alternative-methods

subf="seqs.to.remove" 

rm cluster0000

for f in cluster*
do
	if grep -q 'tig0' $f; then
		cat $f | grep -v '^>' | grep -v 'tig0' |  cut -f 2 | cut -f 2 -d " " | cut -f 2 -d '>' | sed 's/\.\.\.//g' >> $subf
	fi
done

rm cluster*

# fix the seqs to remove file also, remove any trailing "_"
#cat seqs.to.remove | cut -f 1,2,3 -d "_" > temp
#mv temp seqs.to.remove

# get the list of all the sequences for clustering

less $1.fasta | grep '^>' | cut -f 2 -d '>'  > DPTE.names
less $2.fasta  | grep '^>' | cut -f 2 -d '>' > Rmod.names
cat Rmod.names DPTE.names > names.Rmod.DPTE

sort names.Rmod.DPTE > names.Rmod.DPTE.sorted
sort seqs.to.remove > seqs.to.remove.sorted

# Only get the sequences that were not removed using comm -23

comm -23 names.Rmod.DPTE.sorted seqs.to.remove.sorted > seqs.to.keep 


xargs samtools faidx all_consensi.fasta < seqs.to.keep > seqs.tocluster.fasta

/programs/cd-hit-v4.6.1-2012-08-27/cd-hit-est -aS 0.8 -c 0.8 -g 1 -G 0 -A 80 -M 10000 -i seqs.tocluster.fasta -o consensi_clusters2.fasta -T 20

cat consensi_clusters2.fasta $3.fasta > consensi_clusters_final.fasta 




