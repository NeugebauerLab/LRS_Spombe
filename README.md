# LRS_Spombe
This repository includes scripts accompanying the study on RNA porcessing in the fission yeast S. pombe by long-read sequencing

Long-read SMRT cDNA sequencing of nascent RNA from exponentially growing S. pombe 972h- and prp2-1 cells first at 26C
and then for 2 hours at 37C was employed to obtain splicing information from single transcripts. 
Nascent RNA was prepared from the yeast chromatin fraction (Carrillo Oesterreich, Preibisch, Neugebauer, Mol Cell 2010). 
The nascent 3’ end was labeled with a 3’ DNA adaptor through ligation. 
The adaptor sequence served as template for full-length reverse transcription and double-stranded cDNA was obtained 
in a transcriptome-wide PCR. The dataset contains data from 3 replicates per strain. Scripts for demultiplexing the samples, 
filtering for unique barcode matching and removal of transcripts with short polyA tails are deposited here.

Sample barcodes (for dataset published http://dx.doi.org/10.7910/DVN/PW1KEG) are:
sample_id barcode barcode_reverse_complement
SP_wt_BC1	TCAGACGATGCGTCAT	ATGACGCATCGTCTGA
SP_wt_BC2	CTATACATGACTCTGC	GCAGAGTCATGTATAG
SP_wt_BC3	TACTAGAGTAGCACTC	GAGTGCTACTCTAGTA
SP_prp2_BC4	TGTGTATCAGTACATG	CATGTACTGATACACA
SP_prp2_BC5	ACACGCATGACACACT	AGTGTGTCATGCGTGT
SP_prp2_BC6	GATCTCTACTATATGC	GCATATAGTAGAGATC
spike_in_BC7	ACAGTCTATACTGCTG	CAGCAGTATAGACTGT
