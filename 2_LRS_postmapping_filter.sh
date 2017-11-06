
#! /bin/sh

#script to process mapped transcriptome PacBio data, created by Lydia Herzel 08/04/2016
# this script filters for reads which are assigned to multiple samples due to barcode errors
# it retains reads with the same id, but different sequence, originating from concatamers



if [ $# -ne 1 ]
then
	echo "\nThis script filters mapped transcriptome PacBio data for wrong barcode matching."
	echo "Please specify the mapped SAM file names in an extra file, originating from one pooled sample.\n"
	exit  
fi

#input
FILE_NAMES=$1
cat $FILE_NAMES > inSAM_script

cat $(cat $FILE_NAMES | sed 's/\n/ /') > combined.sam

# cut -f 1-4 combined.sam | sort | uniq -c | awk ' $1==1 { print $0 }' | wc -l # 69969
# cut -f 1-4 combined.sam | sort | uniq -c | awk ' $1==2 { print $0 }' | wc -l # 344, 0.5% double assigned

cut -f 1-4 combined.sam | sort | uniq -c | awk ' $1>1 { print $2"\t"$3"\t"$4"\t"$5 }' > non_unique_assignment

# use Rscript to filter for read id + mapping position
# cigar string and seqeunce might differ due to filtering artifacts
Rscript /Users/herzel/Documents/LAB_STUFF/sequencing_array/PacBio_sequencing/SP_Prp2/scripts/LRS_postmapping_filter.R

sed 's/.sam/_filt.sam/' inSAM_script > inSAM_filt

for FILE_NAMES_filt in $(cat inSAM_filt)
	do
		cat $FILE_NAMES_filt | sed 's/^"m/m/' | sed 's/ccs"	/ccs	/' > tmp
		mv tmp $FILE_NAMES_filt
		igvtools index $FILE_NAMES_filt
	done
