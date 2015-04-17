# how to generate the final tables
# chimeric.igv.bam is the results of taking the bwa alignments, marking duplicates and then running the chimeric.py script on them

for file in $(find . -name *.sorted.chimeric.igv.bam); do python ~/cache/li_hiv/scripts/orientation.py $file K03455 > `dirname $file`/`basename $file .sorted.chimeric.igv.bam`.table; done

for file in $(find . -name *.table.bed)
do
	bedtools intersect -wao -a $file -b ~/cache/li_hiv/metadata/features.bed.gz | uniq | cut -f1,2,3,4,5,6,7,11,12,13,14 > `dirname $file`/`basename $file .table.bed`.annotated.bed
done
