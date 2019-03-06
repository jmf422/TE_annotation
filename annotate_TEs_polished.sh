#! /bin/bash

# nohup sh annotate_TEs_polished.sh <species> &

#date
d1=$(date +%s)
echo $1

# do from /workdir/jmf422/Polished_genomes directory

cd $1
mkdir Annotation_files

######################## REPEAT MODELER #########################


mkdir RepeatModeler

cp *.fasta RepeatModeler
cd RepeatModeler

cat *fasta | cut -f 1 -d '-' > $1.assembly.fasta

/programs/RepeatModeler/BuildDatabase -name $1 $1.assembly.fasta

/programs/RepeatModeler/RepeatModeler -pa 12 -database $1

cp RM*/consensi.fa ../Annotation_files/consensi.fasta

echo "RepeatModeler done!"

#date
d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)

cd ..


#################### LTR FINDER #################################

mkdir LTR_finder

cp RepeatModeler/$1.assembly.fasta LTR_finder 

cd /workdir/jmf422/software/LTR_Finder/source

./ltr_finder /workdir/jmf422/Polished_genomes/$1/LTR_finder/$1.assembly.fasta > /workdir/jmf422/Polished_genomes/$1/LTR_finder/$1.LTRfinder.out

echo "LTR finder done!"

# now run LTR retriever to process the output

cd /workdir/jmf422/software/LTR_retriever

./LTR_retriever -genome /workdir/jmf422/Polished_genomes/$1/LTR_finder/$1.assembly.fasta -infinder /workdir/jmf422/Polished_genomes/$1/LTR_finder/$1.LTRfinder.out -threads 12

cd /workdir/jmf422/Polished_genomes/$1/LTR_finder

mv $1.assembly.fasta.LTRlib.fa $1.LTRlib.fasta
cp $1.LTRlib.fasta ../Annotation_files

cd ..
####################### dnaPipeTE ##############################

export JAVA_HOME=/usr/local/jdk1.8.0_121
export PATH=$JAVA_HOME/bin:$PATH

mkdir DPTE
cp RepeatModeler/$1.assembly.fasta DPTE

cd DPTE


art_dir="/workdir/cg629/bin/art_bin_MountRainier"
dnaPipeTE_dir="/workdir/cg629/bin/dnaPipeTE"

cd $art_dir

./art_illumina -ss HS25 -l 150 -f 10 -ef -i /workdir/jmf422/Polished_genomes/$1/DPTE/$1.assembly.fasta -o /workdir/jmf422/Polished_genomes/$1/DPTE/$1.HS25.150.se

cd /workdir/jmf422/Polished_genomes/$1/DPTE

# get the size of the assembly
SIZE=`cat $1.assembly.fasta | grep -v '^>' | wc -m`
CUTOFF=`bc<<<"scale=6; ($SIZE*0.00005)/150"`

cd $dnaPipeTE_dir

python3 dnaPipeTE.py -input /workdir/jmf422/Polished_genomes/$1/DPTE/$1.HS25.150.se.fq -output /workdir/jmf422/Polished_genomes/$1/DPTE/$1.dnaPipeTE_0.25_s2 -cpu 12 -genome_size $SIZE -genome_coverage 0.25 -sample_size 2 -species Drosophila -Trin_glue 2

# now filter the data

cd /workdir/jmf422/Polished_genomes/$1/DPTE
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
