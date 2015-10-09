HIV integration project
======

## background
This is a small set of scripts to identify integration sites and orientation of
viral integrations into a host genome. The strategy is to augment the host
genome with a contig representing the viral genome and look for chimeric
alignments where one part of the alignment is in the host genome and the other
part of the alignment is in the viral genome. We use those alignments to try
to call integration sites. There are some previously existing tools to do this,
such as http://bioinfo.mc.vanderbilt.edu/VirusFinder/ and http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3582262/ and more recently http://www.nature.com/articles/srep11534.
We had some trouble getting dependable results from some of the earlier methods so
we wrote up these scripts which seem to give reasonable results.

## preparing the genome
Augment the human genome with the virus genome. Do that by catting
on the virus genome as an extra chromosome to the human genome,
and then index that augmented genome with bwa.

## aligning a sample
Samples must be demultiplexed. Once they are demultiplexed, align to
the augmented human genome with bwa.

`bwa mem -t $THREADS $GENOME ${sample}_R1_001.fastq.gz ${sample}_R2_001.fastq.gz`

## mark duplicates
mark duplicates with sambamba

`sambamba markdup bwa-mem/$samplename.tmp.bam bwa-mem/$samplename.bam`

## keep only chimeric alignments
This script keeps alignments that are chimeric for the human genome and the
virus genome. In this example we used K03455 as the virus genome. This skips
duplicates and identifies chimeric alignments in two ways. 1) if the alignment
is marked as supplementary and is in the virus sequence and 2) if the alignment
has a supplementary alignment associated with it and one of the alignments
is in the virus sequence and the other is in the human sequence.

`python chimeric.py bwa-mem/$samplename.bam K03455`

## calculate integration and orientation of integrations
python  ~/cache/li_hiv/scripts/orientation.py bam_file virus_contig > out.sites
