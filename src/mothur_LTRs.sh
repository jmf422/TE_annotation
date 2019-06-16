#! /bin/bash
# mothur_LTRs.sh <fasta prefix>

# this script clusters the LTR sequences and chooses the appropriate consensus sequence for each OTU cluster


#date
d1=$(date +%s)
echo $1

# first cluster the LTR sequences - shown to get them to the family level

cat $1.fasta | sed 's/:/xxx/g' > $1.2.fasta 

# make name file

cat $1.2.fasta | grep '^>' | cut -f 2 -d '>' | awk -v OFS="\t" '{print $1,$1}' > $1.names

/programs/mafft/bin/mafft $1.2.fasta > $1.aligned.fasta

/programs/mothur/mothur "#dist.seqs(fasta=$1.aligned.fasta, calc=onegap, countends=T, cutoff=0.20, processors=4)"


# get the value to use
#label=`cat $1.aligned.dist | cut -f 3 | sort -nr | head -n 1`

# now run clustering

/programs/mothur/mothur "#cluster(column=$1.aligned.dist, name=$1.names, method=nearest)"

# need to get the membership of the clusters.

/programs/mothur/mothur "#bin.seqs(list=$1.aligned.nn.list, fasta=$1.2.fasta, name=$1.names, label=0.20)"

#attach the OTU number to the sequence

cat $1.aligned.nn.0.*.fasta | sed 's/\t/-/g' > $1.OTUs.fasta

# needed to shorten things up
cat $1.OTUs.fasta | sed 's/xxx/:/g' > temp.fasta
cat temp.fasta | sed 's/#LTR\/Gypsy//g' | sed 's/#LTR\/Copia//g' | sed 's/#LTR\/unknown//g' > $1.OTUs.fasta

# now run Rmasker - do not consider the candidates that have similarity to LINE or DNA element

/programs/RepeatMasker/RepeatMasker -lib zebrep.fa -nolow -div 20 -pa 8 $1.OTUs.fasta # changed this

#get the singletons - they will automatically be the consensus

cat $1.OTUs.fasta | grep '^>' | cut -f 2 -d '-' | uniq -c | sed -e 's/^[ \t]*//' | tr -s " " | sed 's| |\t|g' | awk '$1==1 {print $2}' > singleton.OTUs
# grep the group
cat $1.OTUs.fasta | grep '^>' | cut -f 2 -d '>' > $1.OTUs.names
grep -f singleton.OTUs $1.OTUs.names > singleton_seqs

# this might need to be fixed - DNA, LINE
# need to first remove the LTRs that have some extraneous sequence attached.
less $1.OTUs.fasta.out | sed -e 1,3d | sed -e 's/^[ \t]*//' | tr -s " " | sed 's| |\t|g' | cut -f 5,11 > $1.OTUs.fasta.classification

cat $1.OTUs.fasta.classification | grep -E "(DNA|hAT|Mariner)" | cut -f 1 > LTRs.with.insertions

# only want a fasta file with the LTRs that are not in the LTRs.with.insertions file

sort LTRs.with.insertions > LTRs.with.insertions.sorted
sort $1.OTUs.names > LTRs.all.sorted
comm -23 LTRs.all.sorted  LTRs.with.insertions.sorted > LTRs.candidates 

sort singleton_seqs > singleton_seqs.sorted
sort LTRs.candidates > LTRs.candidates.sorted

### new part ####

comm -23 LTRs.candidates.sorted singleton_seqs.sorted > refiner.candidates

xargs samtools faidx $1.OTUs.fasta < refiner.candidates > refiner.candidates.fasta

# get the singleton sequences to combine with

xargs samtools faidx $1.OTUs.fasta < singleton_seqs.sorted > singletons.fasta

# go through each OTU that has multiple members and run refiner

otus=`cat refiner.candidates | cut -f 2 -d '-' | sort -u`

subf="rep.sequences"
for o in $otus
do
	cat refiner.candidates.fasta | grep $o | cut -f 2 -d '>' > $o.seqs
	type=`head -n 1 $o.seqs | cut -f 2 -d '_' | cut -f 1 -d '-'`
	xargs samtools faidx refiner.candidates.fasta < $o.seqs > $o.seqs.fasta
	perl /workdir/jmf422/software/RepeatModeler_new/dist/Refiner $o.seqs.fasta
	rm $o.seqs
	rm $o.seqs.fasta
	sed "s/>.*/>chr_refiner_$type-$o/" $o.seqs.fasta.refiner_cons > $o.refiner.cons.fasta
	rm $o.seqs.fasta.refiner_cons
done 

cat *.refiner.cons.fasta singletons.fasta > clustered_refined_LTRlib.fasta


#date
d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)