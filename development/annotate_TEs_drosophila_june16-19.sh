#! /bin/bash

# nohup sh annotate_TEs_drosophila_june16-19.sh <species> &

#date
d1=$(date +%s)
echo $1

# do from /workdir/jmf422/Final


mkdir Annotation_files

######################## REPEAT MODELER #########################


mkdir RepeatModeler

cat *.fasta | cut -f 1 -d ' ' > $1.assembly.fasta

cp $1.assembly.fasta RepeatModeler
cd RepeatModeler


/workdir/jmf422/software/RepeatModeler_new/dist/BuildDatabase -name $1 $1.assembly.fasta

/workdir/jmf422/software/RepeatModeler_new/dist/RepeatModeler -pa 16 -database $1

cp RM*/consensi.fa ../Annotation_files

echo "RepeatModeler done!"

#date
d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)

cd ..

##################### LTR Harvest #############################

mkdir LTR_Harvest

cp RepeatModeler/$1.assembly.fasta LTR_Harvest

cd LTR_Harvest

/workdir/jmf422/software/genometools-1.5.9/bin/gt suffixerator -db $1.assembly.fasta -indexname $1.assembly -tis -suf -lcp -des -ssp -sds -dna

/workdir/jmf422/software/genometools-1.5.9/bin/gt ltrharvest -index $1.assembly > $1.ltrharvest.out

echo "LTR harvest done!"

# now run LTR retriever to process the output


/workdir/jmf422/software/LTR_retriever/LTR_retriever -genome $1.assembly.fasta -inharvest $1.ltrharvest.out -threads 12


mv $1.assembly.fasta.LTRlib.fa $1.LTRlib.fasta
cp $1.LTRlib.fasta ../Annotation_files

#date
d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)

cd ..

####################### dnaPipeTE ##############################

export JAVA_HOME=/usr/local/jdk1.8.0_121
export PATH=$JAVA_HOME/bin:$PATH

mkdir DPTE
cp RepeatModeler/$1.assembly.fasta DPTE

cd DPTE


art_dir="/workdir/cg629/bin/art_bin_MountRainier"
dnaPipeTE_dir="/workdir/cg629/bin/dnaPipeTE"


$art_dir/art_illumina -ss HS25 -l 150 -f 10 -ef -i $1.assembly.fasta -o $1.HS25.150.se


# get the size of the assembly
SIZE=`cat $1.assembly.fasta | grep -v '^>' | wc -m`
CUTOFF=`bc<<<"scale=6; ($SIZE*0.00005)/150"`


python3 $dnaPipeTE_dir/dnaPipeTE.py -input $1.HS25.150.se.fq -output $1.dnaPipeTE_0.25_s2 -cpu 12 -genome_size $SIZE -genome_coverage 0.25 -sample_size 2 -species Drosophila -Trin_glue 2 # only important for classification?

# now filter the data

cd /workdir/jmf422/$1/DPTE
cp $1.dnaPipeTE_0.25_s2/reads_per_component_and_annotation .
cp $1.dnaPipeTE_0.25_s2/Trinity.fasta .

cat reads_per_component_and_annotation | awk -v CUTOFF="$CUTOFF" '$1>CUTOFF {print $0}' | cut -f 3 -d " " > $1.DPTE.contigs

# now get these sequences

xargs samtools faidx Trinity.fasta < $1.DPTE.contigs > $1.DPTE.contigs.fasta

cp $1.DPTE.contigs.fasta ../Annotation_files

echo "done dnaPipeTE!"

#date
d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)