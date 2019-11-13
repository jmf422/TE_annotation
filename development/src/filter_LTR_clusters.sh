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
	if grep -q 'utg' $f; then
		cat $f | grep -v '^>' | grep -v 'chr' |  cut -f 2 | cut -f 2 -d " " | cut -f 2 -d '>' | sed 's/\.\.\.//g' >> $subf
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


csplit --digits=4 --quiet --prefix=clstr consensi_clusters2.fasta.clstr "/^>/" "{*}"

rm clstr0000

# go through the clstr file 
# get the ones that are not singletons and run refiner

for f in clstr*
do
	if grep -q '^1' $f; then
		cat $f | head -n 1 | cut -f 2 -d '>' >> not.singletons
		cat $f | sed '1d' | cut -f 2 -d '>'  | cut -f 1 -d '.' > $f.candidates
		xargs samtools faidx seqs.tocluster.fasta < $f.candidates > $f.candidates.fasta
		perl /workdir/jmf422/software/RepeatModeler_new/dist/Refiner $f.candidates.fasta
		rm $f.candidates
		rm $f
		sed "s/>.*/>fam_refiner_$f/" $f.candidates.fasta.refiner_cons > $f.refiner.cons.fasta
		rm $f.candidates.fasta.refiner_cons
	fi
done

# get the sequences of the singletons now

for f in clstr*
do
	cat $f | sed '1d' | cut -f 2 -d '>'  | cut -f 1 -d '.' >> singleton.seqs
done


xargs samtools faidx seqs.tocluster.fasta < singleton.seqs > singleton.seqs.fasta

cat *refiner.cons.fasta singleton.seqs.fasta > nonLTRs.cdhit.refiner.fasta

cat nonLTRs.cdhit.refiner.fasta clustered_refined_LTRlib.fasta > consensi_clusters_final.fasta

rm clstr*

#date
d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)



