#! /bin/Rscript

#Rscript to eliminate transcripts, which originate from cleaved transcripts with short polyA tails

## INPUT
sam <- read.table("input_short", header = F, stringsAsFactors = F)
#sam <- read.table("/Users/herzel/Documents/LAB_STUFF/sequencing_array/PacBio_sequencing/compilation_for_publication/mapped_data/polyA_filt/input_short", header = F, stringsAsFactors = F)

names(sam) <- c("strand", "chromosome", "start", "cigar", "sequence")
#sam <- sam[sample(nrow(sam))[1:10000],]
bed <- read.table("/Users/herzel/Documents/LAB_STUFF/annotations/Sp_annotations/Sp_ASM294v2_31_Eser_combined_100us_100dsPolyA.bed",header=F,stringsAsFactors=F)

#plus strand in bed + and in sam 0
bed <- split(bed,bed[,1])

#identify read chromosome, search in according bed file for overlapping TSS
#do a similar procedure for gene end
#here read end is searched..

entry_end <- c(1:nrow(sam))
for (n in 1:nrow(sam)) {
  i <- sam[n,]
  
  chr <- i$chromosome
  if(i$strand==0) { strand <- "+" } else { strand <- "-" }
  if(i$strand==16) { start <- as.numeric(i$start) 
  
  } else { 
    
    if (grepl("D",i$cigar)==F & grepl("N",i$cigar)==T) {
      
      N <- strsplit(i$cigar,"N")[[1]]
      fieldN <- grepl("[1-9]$",N)
      N <- sum(as.numeric(sub(".+[A-Z]","",N[fieldN])))
      
      start <- as.numeric(i$start)+nchar(i$sequence)+N
      
    } else if (grepl("D",i$cigar)==T & grepl("N",i$cigar)==T ) { 
      
      D <- strsplit(i$cigar,"D")[[1]]
      fieldD <- grepl("[1-9]$",D)
      D <- sum(as.numeric(sub(".+[A-Z]","",D[fieldD])))
      
      N <- strsplit(i$cigar,"N")[[1]]
      fieldN <- grepl("[1-9]$",N)
      N <- sum(as.numeric(sub(".+[A-Z]","",N[fieldN])))
      
      start <- as.numeric(i$start)+nchar(i$sequence)-D+N
      
    } else if (grepl("D",i$cigar)==T & grepl("N",i$cigar)==F ) { 
      D <- strsplit(i$cigar,"D")[[1]]
      fieldD <- grepl("[1-9]$",D)
      D <- sum(as.numeric(sub(".+[A-Z]","",D[fieldD])))
      
      start <- as.numeric(i$start)+nchar(i$sequence)-D
      
    } else if (grepl("D",i$cigar)==F & grepl("N",i$cigar)==F ) { 
      
      start <- as.numeric(i$start)+nchar(i$sequence)
      
    }
  }
  
  chrRegion <- bed[names(bed)==chr][[1]]
  chrRegion <- chrRegion[chrRegion[,6]==strand,c(2:3)]
  hit <- length(chrRegion[start >=  chrRegion[,1] & start <=  chrRegion[,2],1])
  entry_end[n] <- hit
}

entry_end[entry_end==0] <- F
entry_end[entry_end!=0] <- T

#check, if read end contains stretch of As
polyA <- sam$strand==0 & grepl("[ACGT]AAAA$",sam$sequence) | grepl("A[ACGT]AAA$",sam$sequence) | grepl("AA[ACGT]AA$",sam$sequence) | grepl("AAA[ACGT]A$",sam$sequence) | grepl("AAAA[ACGT]$",sam$sequence) 
polyT <- sam$strand==16 & grepl("^[ACGT]TTTT",sam$sequence) | grepl("^T[ACGT]TTT",sam$sequence) | grepl("^TT[ACGT]TT",sam$sequence) | grepl("^TTT[ACGT]T",sam$sequence) | grepl("^TTTT[ACGT]",sam$sequence)
truncA <- sam$strand==0 & grepl("S$",sam$cigar) & grepl("A$",sam$sequence) 
truncT <- sam$strand==16 & grepl("^[0-9]*S",sam$cigar) & grepl("^T",sam$sequence) 

#overlap entry_end and polyA should be filtered out
filterOut <- entry_end==1 & (polyA==T | polyT==T | truncA==T | truncT==T)
filterOut[filterOut==T] <- 1
write.table(filterOut, "polyA_filt", col.names=F, row.names=F, quote=F, sep="\t")
#write.table(filterOut, "/Users/herzel/Documents/LAB_STUFF/sequencing_array/PacBio_sequencing/compilation_for_publication/mapped_data/polyA_filt/polyA_filt", col.names=F, row.names=F, quote=F, sep="\t")
