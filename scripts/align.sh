#!/bin/bash
GENOME=/net/hsphfs1/srv/export/hsphfs1/share_root/chb/projects/li_hiv/virusscan/genomes/hg19_hiv/seq/hg19_hiv.fa
THREADS=8
DATADIR=/net/hsphfs1/srv/export/hsphfs1/share_root/chb/projects/li_hiv/data/09-2014
mkdir -p bwa-mem
for sample in `python /net/hsphfs1/srv/export/hsphfs1/share_root/chb/projects/li_hiv/code/li_hiv/scripts/get_samples.py $DATADIR`
do
    echo "Submitting $sample"
    samplename=`basename $sample`
    sbatch -J $sample --cpus-per-task=$THREADS --time=1:00:00 --mem=10000 --wrap="bwa mem -t $THREADS $GENOME ${sample}_R1_001.fastq ${sample}_R2_001.fastq | samtools view -Su - | samtools sort - bwa-mem/$samplename.tmp; sambamba markdup bwa-mem/$samplename.tmp.bam bwa-mem/$samplename.bam"
done
