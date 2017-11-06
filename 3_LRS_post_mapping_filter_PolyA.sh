#script to process transcriptome PacBio data

#! /bin/sh

#script to process transcriptome PacBio data, created by Lydia Herzel 03/14/2015
if [ $# -ne 2 ]
then
	echo "\nThis script filters S. pombe PacBio data with respect to the removal of residual polyadenylated transcripts. After it converts the sam-file to a bam & bed-file for splicing analysis"
	echo "Please specify the list of SAM filepaths and the path to the bam-header.\n"
	exit  
fi

#input
FILES=$1
HEADER=$2

for SAM in $(cat $FILES)
	do
		NAME=$(echo $SAM |  awk ' BEGIN {FS="/"}; {print $NF}' | sed 's/.sam//' )
		echo $NAME
		cut -f 2,3,4,6,10 $SAM > input_short
		
		echo "identifying reads with short polyA tails and 3' ends close to annotated polyA sites..."
		Rscript ~/Documents/LAB_STUFF/sequencing_array/PacBio_sequencing/compilation_for_publication/scripts/filter_for_short_polyAtails.R 
		echo "done"
		echo "removing classified reads with A-tail and converting to bam-file"
		paste $SAM polyA_filt | grep "	0$" | cut -f 1-21 > tmp.sam
		cat $HEADER tmp.sam > "$NAME"_filtPolyA.sam
		samtools view -S -b "$NAME"_filtPolyA.sam > tmp.bam
		samtools sort tmp.bam > "$NAME"_filtPolyA.bam
		samtools index "$NAME"_filtPolyA.bam
		bamToBed -bed12 -i "$NAME"_filtPolyA.bam > "$NAME"_filtPolyA.bed
		echo "done"
	done
	
rm tmp.sam tmp.bam input_short polyA_filt