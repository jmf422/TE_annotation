#! /bin/bash
# this file evaluates which consensi mask the genome with RepeatMasker

# sh eval_consensi_Rmasker.sh <library> <genome> <species>

#date
d1=$(date +%s)
echo $1

mkdir mask_eval

cp $1.fasta mask_eval
cd mask_eval

cp /workdir/jmf422/zebrafish/RepeatModeler/$2.fasta .

# first run Repeat masker

/programs/RepeatMasker/RepeatMasker -lib $1.fasta -nolow -pa 8 $2.fasta 

# then process the file

less $2.fasta.out | sed -e 1,3d | sed -e 's/^[ \t]*//' | tr -s " " | sed 's| |\t|g' | cut -f 5,6,7,10-14 > $2.libmasked.bed

cat $2.libmasked.bed | cut -f 4 | sort | uniq -c | sed -e 's/^[ \t]*//' | sed 's/ /\t/g' > $2.libmasked.uniq

cat $2.libmasked.uniq | awk '$1>=2 {print $2}' > $2.libmasked.consensi

cat $1.fasta | grep '^>' | cut -f 2 -d ">" | sort > $1.sorted

comm -23 $1.sorted $2.libmasked.consensi > $1.not.masked # these are the consensi not masked 
# remove them from the list

comm -23 $1.sorted $1.not.masked > $1.filtered

# get fasta file now

xargs samtools faidx $1.fasta < $1.filtered > $1.filtered.fasta

cp  $1.filtered.fasta ..


d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)