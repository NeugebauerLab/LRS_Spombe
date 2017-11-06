
#! /bin/sh

# script to process transcriptome PacBio data, created by Lydia Herzel (herzel@mit.edu, 07/29/2017)
# this script filters for sample barcodes and deconcatenates reads with multiple cDNA inserts

# the reads get filtered for a sample barcode at the 5' and 3' end, then trimmed for the adaptor, filtered for duplicates,
# trimmed for the 5N random barcode
# reverse complements sequence
# maps data

if [ $# -ne 3 ]
then
	echo "\nThis script executes processes transcriptome PacBio data."
	echo "Please specify the FASTQ filepath,the preferred GENOME (Scer/Spombe_EF2) and the file containing two fields (id\tbarcode-sequence).\n"
	exit  
fi

#input
FASTQ=$1
GENOME=$2
SAMPLE_BARCODE=$3
ERROR_RATE1=0.1
ERROR_RATE2=0.12

#filter for sample barcode
# 1st filter for 3adaptor 

while read lineA           
	do
		BC=$(echo $lineA | awk ' { print $2 } ')
		revBC=$(echo $BC | rev | sed 's/A/1/g' | sed 's/C/2/g' | sed 's/G/3/g' | sed 's/T/4/g' | sed 's/1/T/g' | sed 's/2/G/g' | sed 's/3/C/g' | sed 's/4/A/g' )
		ID=$(echo $lineA | awk ' { print $1 } ')

		echo $BC"\t"$revBC

		cutadapt -a CTGTAGGCACCATCAATGTACTCTGCGTTGATACCACTGCTT"$revBC" -n 1 -O 58 -e $ERROR_RATE1 -m 10 --match-read-wildcards --untrimmed-output=untrimmed.fastq --info-file=info1f.fastq  -r trimmed1f_rest.fastq -o trimmed1f.fastq $FASTQ

		cat untrimmed.fastq | tr '\n' '	' | sed 's/@m1/\
@m1/g' | sed '/^$/d' > untrimmed.tab

		cut -f 1 untrimmed.tab > ID
		cut -f 2 untrimmed.tab | rev | sed 's/A/1/g' | sed 's/C/2/g' | sed 's/G/3/g' | sed 's/T/4/g' | sed 's/1/T/g' | sed 's/2/G/g' | sed 's/3/C/g' | sed 's/4/A/g' > SEQ
		cut -f 3 untrimmed.tab > PLUS
		cut -f 4 untrimmed.tab | rev > QUAL

		paste ID SEQ PLUS QUAL | sed 's/	/\
/g' > untrimmed_rev.fastq

		rm ID SEQ PLUS QUAL untrimmed.tab

		cutadapt -a CTGTAGGCACCATCAATGTACTCTGCGTTGATACCACTGCTT"$revBC" -n 1 -O 58 -e $ERROR_RATE1 -m 10 --match-read-wildcards --untrimmed-output=untrimmed2.fastq --info-file=info2f.fastq  -r trimmed2f_rest.fastq -o trimmed2f.fastq untrimmed_rev.fastq
		# untrimmed2.fasta mainly rRNA, incomplete adaptor etc #4506, 9%

		# generate fastQ of rest file, question: where are adaptors trimmed, most 5' and most 3' or both most 5'?
		# reverse adaptor etc might be still present as concatenation is not strand specific.. indeed n=3613 perfect matches, fwd adaptor 400 matches!
		# proceed with sample_bc trimming, proper orientation of barcode & SMARTer adaptor should go with previous trimming
		# again rest is important to keep for potential other cDNAs

		cat trimmed1f.fastq trimmed2f.fastq > trimmedf1.fastq
		rm trimmed1f.fastq trimmed2f.fastq untrimmed_rev.fastq
		cat info1f.fastq info2f.fastq > infof1.fastq
		rm info1f.fastq info2f.fastq 
		cat trimmed1f_rest.fastq trimmed2f_rest.fastq > trimmedf1_rest.fastq
		rm trimmed1f_rest.fastq trimmed2f_rest.fastq 

# do second round of filtering for potential concatamers, which have not been isolated as reads
		cutadapt -a CTGTAGGCACCATCAATGTACTCTGCGTTGATACCACTGCTT"$revBC" -n 1 -O 58 -e $ERROR_RATE2 -m 10 --match-read-wildcards --untrimmed-output=untrimmed3.fastq --info-file=info3f.fastq  -r trimmed3f_rest.fastq -o trimmed3f.fastq untrimmed2.fastq

		cat untrimmed3.fastq | tr '\n' '	' | sed 's/@m1/\
@m1/g' | sed '/^$/d' > untrimmed.tab

		cut -f 1 untrimmed.tab > ID
		cut -f 2 untrimmed.tab | rev | sed 's/A/1/g' | sed 's/C/2/g' | sed 's/G/3/g' | sed 's/T/4/g' | sed 's/1/T/g' | sed 's/2/G/g' | sed 's/3/C/g' | sed 's/4/A/g' > SEQ
		cut -f 3 untrimmed.tab > PLUS
		cut -f 4 untrimmed.tab | rev > QUAL

		paste ID SEQ PLUS QUAL | sed 's/	/\
/g' > untrimmed_rev.fastq

		rm ID SEQ PLUS QUAL untrimmed.tab

		cutadapt -a CTGTAGGCACCATCAATGTACTCTGCGTTGATACCACTGCTT"$revBC" -n 1 -O 58 -e $ERROR_RATE2 -m 10 --match-read-wildcards --untrimmed-output=untrimmed4.fastq --info-file=info4f.fastq  -r trimmed4f_rest.fastq -o trimmed4f.fastq untrimmed_rev.fastq
		cat trimmed3f.fastq trimmed4f.fastq > trimmedf2.fastq
		rm trimmed3f.fastq trimmed4f.fastq untrimmed_rev.fastq
		cat info3f.fastq info4f.fastq > infof2.fastq
		rm info3f.fastq info3f.fastq 
		cat trimmed3f_rest.fastq trimmed4f_rest.fastq > trimmedf2_rest.fastq
		rm trimmed3f_rest.fastq trimmed4f_rest.fastq 

# bring both rounds of trimming together
		cat trimmedf1.fastq trimmedf2.fastq > trimmedf.fastq
		
		# trim for 5' end barcode - SMARTer sequence & NNGGG from template switching (modification inserted 07/29/2016)
#		cutadapt -g "$BC"AAGCAGTGGTATCAACGCAGAGTACNNGGG -n 1 -O 46 -e $ERROR_RATE1 -m 10 --match-read-wildcards --discard-untrimmed  --info-file=info2f_"$ID".fastq  -r trimmed2f_rest_"$ID".fastq -o trimmed2f_"$ID".fastq trimmedf.fastq
		cutadapt -g "$BC"AAGCAGTGGTATCAACGCAGAGTACNNNNN -n 1 -O 46 -e $ERROR_RATE1 -m 10 --match-read-wildcards --discard-untrimmed  --info-file=info2f_"$ID".fastq  -r trimmed2f_rest_"$ID".fastq -o trimmed2f_"$ID".fastq trimmedf.fastq

		echo "Rscript filterFASTQforPCRduplicates launched .."
		cp trimmed2f_"$ID".fastq input.fastq
		Rscript /Users/herzel/Documents/LAB_STUFF/sequencing_array/PacBio_sequencing/scripts/filterFASTQforPCRduplicates_3prime5N.R
		echo ".. Rscript filterFASTQforPCRduplicates finished"
		
		#trim by 5N at 3'end
		[ -f output.fastq ] && fastx_trimmer -Q 33 -t 5 -i output.fastq -o trimmed2f_noDupl_"$ID".fastq || echo "something's wrong .. output.fastq not found"
		rm output.fastq
		
		fastx_trimmer -Q 33 -t 5 -i input.fastq -o trimmed2f_noDupl_"$ID".fastq
		
		echo $GENOME
		if [ "$GENOME" == "Spombe_EF2" ]
			then
				echo "Spombe mapping..."
#				gmap -d Spombe_EF2 --min-intronlength=30 --intronlength=850 --localsplicedist=850 --totallength=850 --trimendexons=0 --microexon-spliceprob=0.5 --direction=auto --find-shifted-canonical --allow-close-indels=2 --npaths=1 --nofails --fails-as-input --mapboth -A --format=samse trimmed2f_noDupl_"$ID".fastq > "$ID".sam
				/Users/herzel/programs/gmap-2016-07-11/src/gmap -D /Users/herzel/programs/gmap-2016-07-11/database -d Spombe_EF2 --min-intronlength=20 --max-intronlength-middle=1000 --max-intronlength-ends=1000 --trim-end-exons=0 --localsplicedist=1000 --no-chimeras --allow-close-indels=2 --microexon-spliceprob=0.5 --npaths=1 --nofails --failed-input=no_mapping --mapboth -A --format=samse trimmed2f_noDupl_"$ID".fastq  > "$ID".sam

				echo "finished mapping"
		elif [ "$GENOME" == "Scer" ]
			then
				echo "Scer3 mapping..."
#				gmap -d Scer --min-intronlength=30 --intronlength=1010 --localsplicedist=1010 --totallength=1010 --trimendexons=0 --microexon-spliceprob=0.5 --direction=auto --find-shifted-canonical --allow-close-indels=2 --npaths=1 --nofails --fails-as-input --mapboth -A --format=samse trimmed2f_noDupl_"$ID".fastq  > "$ID".sam
				/Users/herzel/programs/gmap-2016-07-11/src/gmap -D /Users/herzel/programs/gmap-2016-07-11/database -d Scer3 --min-intronlength=20 --max-intronlength-middle=1000 --max-intronlength-ends=1000 --trim-end-exons=0 --localsplicedist=1000 --no-chimeras --allow-close-indels=2 --microexon-spliceprob=0.5 --npaths=1 --nofails --failed-input=no_mapping --mapboth -A --format=samse trimmed2f_noDupl_"$ID".fastq  > "$ID".sam
				echo "finished mapping"
		elif [ "$GENOME" == "PDR5" ]
			then
				echo "PDR5 spike-in mapping..."
		#		gmap -d PDR5 --direction=auto --find-shifted-canonical --allow-close-indels=2 --npaths=1 --nofails --fails-as-input --mapboth -A --format=samse trimmed2f_noDupl_"$ID".fastq  > Spike_in_"$ID".sam
				# currently soft clipping... does not reflect correct 5'end
		#		gsnap -d PDR5 --trim-mismatch-score=0 trimmed2f_noDupl_"$ID".fastq  > Spike_in_"$ID".sam # gsanp takes only fasta
				/Users/herzel/programs/gmap-2016-07-11/src/gmap -D /Users/herzel/programs/gmap-2016-07-11/database -d PDR5 --min-intronlength=20 --max-intronlength-middle=1000 --max-intronlength-ends=1000 --trim-end-exons=0 --localsplicedist=1000 --no-chimeras --allow-close-indels=2 --microexon-spliceprob=0.5 --npaths=1 --nofails --failed-input=no_mapping --mapboth -A --format=samse trimmed2f_noDupl_"$ID".fastq  > "$ID".sam
				
				echo "finished mapping"
				
				grep -v "^>" Spike_in_"$ID".sam | grep -v "^[ACGT]" > tmp.sam
				igvtools sort tmp.sam Spike_in_"$ID".sorted.sam
				cat Spike_in_"$ID".sorted.sam | awk ' $5==40 { print $0} ' > Spike_in_"$ID".sam
				igvtools index Spike_in_"$ID".sam
				rm Spike_in_"$ID".sorted.sam		
		fi

		grep -v "^>" "$ID".sam | grep -v "^[ACGT]" > tmp.sam
		igvtools sort tmp.sam "$ID".sorted.sam
		cat "$ID".sorted.sam | awk ' $5==40 { print $0} ' > "$ID".sam
		igvtools index "$ID".sam
		rm "$ID".sorted.sam		
		
	done<$SAMPLE_BARCODE