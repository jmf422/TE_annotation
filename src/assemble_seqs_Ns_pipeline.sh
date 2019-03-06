#! /bin/bash

#assemble_seqs_Ns_pipeline.sh <Ns file> <genome file> <species>
# RmodConsensi.Ns dmel.MHAP.all

# date
d1=$(date +%s)
echo $HOSTNAME

cp /workdir/jmf422/Polished_genomes/$3/RepeatModeler/$2.fasta .
/programs/RepeatMasker/RepeatMasker -lib $1.fasta -nolow -pa 12 $2.fasta 

mv $2.fasta.out Genome.Ns.masking.out

cat Genome.Ns.masking.out | sed -e 1,3d | sed -e 's/^[ \t]*//' | tr -s " " | sed 's| |\t|g' | awk -v OFS="\t" '{print $5,$6,$7,$10,$11,$12,$13,$14,$2}' | awk '$9<10 {print $0}'  > Rmasker.bed


# get each unique sequence name 
cat Rmasker.bed | cut -f 4 | sort -u > seqs.masked

# loop through each one and assemble


while read line
do
	cat Rmasker.bed | grep -w $line | cut -f 1,2,3 > $line.bed
	bedtools getfasta -fi $2.fasta -bed $line.bed -fo $line.fasta
	/programs/CAP3/cap3i $line.fasta 
	rm $line.bed
	csplit --digits=4 --quiet --prefix=contig_num $line.fasta.cap.contigs.qual "/>/" "{*}"
	rm contig_num0000
	for f in contig_num*
	do
		contig=`head -n 1 $f | cut -f 2 -d '>'`
		avg_qual=`cat $f | sed '1d' | sed 's/ /\n/g' |  sed '/^$/d' | awk '{sum += $1} END {print sum/NR}'`
		length=`cat $f | sed '1d' | sed 's/ /\n/g' |  sed '/^$/d' | wc -l`
		printf "%s\t%f\t%i\n" $contig $avg_qual $length >> $line.contigs
	done
	cat $line.contigs | sort -nrk2,2 -nrk3,3 | head -n 1 | cut -f 1 > $line.best.contig
	xargs samtools faidx $line.fasta.cap.contigs < $line.best.contig > $line.best.contig.fasta
	cat $line.best.contig.fasta | sed "s/Contig.*/${line}/g" > $line.best.contig2.fasta
	rm contig_num0*
done < seqs.masked	

cat *best.contig2.fasta > Nseq.redo.contigs.fasta

# choose the contigs with high quality next

while read line
do
	qual=`cat $line.contigs | sort -nrk2,2 -nrk3,3 | head -n 1 | cut -f 2`
	len=`cat $line.contigs | sort -nrk2,2 -nrk3,3 | head -n 1 | cut -f 3`
	printf "%s\t%f\t%i\n" $line $qual $len >> contig.qual.lengths
done < seqs.masked

cat contig.qual.lengths | awk '{if($2>=50 && $3>200 && $3<10000) print $0}' | cut -f 1 > contigs.passed
xargs samtools faidx Nseq.redo.contigs.fasta < contigs.passed > Nseq.contigs.passed.fasta 

mkdir Ns_files
mv rnd* Ns_files



#date
d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)