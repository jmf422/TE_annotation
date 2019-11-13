#! /bin/bash
# sh filter_long_only.sh <fasta file prefix>  


#date
d1=$(date +%s)
echo $1

# first, calculate the length of the sequences

samtools faidx $1.fasta
cat $1.fasta.fai | cut -f 1,2 > $1.lengths

cat $1.lengths | awk '$2<=10000 {print $1}' > $1.goodlength
cat $1.lengths | awk '$2>10000 {print $1}' > $1.long

# now get the corresponding fasta sequences

xargs samtools faidx $1.fasta < $1.goodlength > $1.goodlength.fasta

xargs samtools faidx $1.fasta < $1.long > $1.long.fasta


#date
d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)