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

# now run clustering

/programs/mothur/mothur "#cluster(column=$1.aligned.dist, name=$1.names, method=nearest)"

# need to get the membership of the clusters.

/programs/mothur/mothur "#bin.seqs(list=$1.aligned.nn.list, fasta=$1.2.fasta, name=$1.names, label=0.20)"

#attach the OTU number to the sequence

cat $1.aligned.nn.0.20.fasta | sed 's/\t/-/g' > $1.OTUs.fasta

# needed to shorten things up
cat $1.OTUs.fasta | sed 's/xxx/:/g' > temp.fasta
cat temp.fasta | sed 's/#LTR\/Gypsy//g' | sed 's/#LTR\/Copia//g' | sed 's/#LTR\/unknown//g' > $1.OTUs.fasta

# now run Rmasker - do not consider the candidates that have similarity to LINE or DNA element

/programs/RepeatMasker/RepeatMasker -species Drosophila -nolow -div 20 $1.OTUs.fasta

#get the singletons - they will automatically be the consensus

cat $1.OTUs.fasta | grep '^>' | cut -f 2 -d '-' | uniq -c | sed -e 's/^[ \t]*//' | tr -s " " | sed 's| |\t|g' | awk '$1==1 {print $2}' > singleton.OTUs
# grep the group
cat $1.OTUs.fasta | grep '^>' | cut -f 2 -d '>' > $1.OTUs.names
grep -f singleton.OTUs $1.OTUs.names > singleton_seqs

# need to first remove the LTRs that have some extraneous sequence attached.
less $1.OTUs.fasta.out | sed -e 1,3d | sed -e 's/^[ \t]*//' | tr -s " " | sed 's| |\t|g' | cut -f 5,11 > $1.OTUs.fasta.classification

cat $1.OTUs.fasta.classification | grep -E "(DNA|LINE)" | cut -f 1 > LTRs.with.insertions

# only want a fasta file with the LTRs that are not in the LTRs.with.insertions file

sort LTRs.with.insertions > LTRs.with.insertions.sorted
sort $1.OTUs.names > LTRs.all.sorted
comm -23 LTRs.all.sorted  LTRs.with.insertions.sorted > LTRs.candidates 


# combine with the singletons
cat LTRs.candidates singleton_seqs | sort -u > LTRs.candidates.all

# now get fasta file

xargs samtools faidx $1.OTUs.fasta < LTRs.candidates.all > $1.OTUs.candidates.fasta

# now get the lengths of the file

samtools faidx $1.OTUs.candidates.fasta
cut -f 1,2 $1.OTUs.candidates.fasta.fai > $1.OTUs.candidates.lengths

otus=`cat $1.OTUs.candidates.lengths | cut -f 1 | cut -f 2 -d '-' | sort -u`

# choose the longest sequence as the OTU representatitive
subf="rep.sequences"
for o in $otus
do
	cat $1.OTUs.candidates.lengths | grep $o | sort -nrk2,2 | head -n 1 | cut -f 1 >> $subf
done 

# now get the fasta sequence of the rep 

xargs samtools faidx $1.OTUs.fasta < rep.sequences > $1.LTR.OTUreps.fasta

mv *lengths *candidates* LTRs* *names *singleton*  extra_files

#date
d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)