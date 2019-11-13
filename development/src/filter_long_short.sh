#! /bin/bash
# sh filter_long_short.sh <fasta file prefix> <short threshold>


#date
d1=$(date +%s)
echo $1

CUTOFF=$2
echo $CUTOFF
# first, calculate the length of the sequences

samtools faidx $1.fasta
cat $1.fasta.fai | cut -f 1,2 > $1.lengths

cat $1.lengths| awk '$2<=10000 {print $0}' | awk -v CUTOFF="$CUTOFF" '$2>=CUTOFF {print $1}' > $1.goodlength

cat $1.lengths | awk '$2>10000 {print $1}' > $1.long
cat $1.lengths | awk -v CUTOFF="$CUTOFF" '$2<CUTOFF {print $1}' > $1.short
cat $1.long $1.short > $1.badlength

# now get the corresponding fasta sequences

xargs samtools faidx $1.fasta < $1.goodlength > $1.goodlength.fasta

xargs samtools faidx $1.fasta < $1.badlength > $1.badlength.fasta


#date
d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)