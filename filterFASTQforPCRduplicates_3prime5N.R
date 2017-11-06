#!/bin/Rscript

library(ShortRead)

#read in FastQ file (already processed with cutadapt)
readFastq("input.fastq") -> fastQ

#filter for unique sequences: 
#sort reads
as.character(as.data.frame(sread(fastQ))[,1]) -> reads
order(reads) -> nrow
reads <- reads[nrow]

reads_length <- nchar(reads)
nachfolger <- c(reads_length[2:length(reads_length)],reads_length[length(reads_length)])

window <- 15

#lapply(,rev) to get end of read
reads_10 <- lapply(lapply(strsplit(reads,""),rev),'[',1:window)
nachfolger_10 <- c(reads_10[2:length(reads_10)],reads_10[length(reads_10)])

reads_10 <- mapply(reads_10, nachfolger_10, FUN=list, SIMPLIFY=FALSE)
non_uniq_10 <- lapply(lapply(reads_10,function(i) {i[[1]]==i[[2]]}),table)
non_uniq_10 <- sapply(non_uniq_10,function(i) {i[names(i)==T]==window})

#alternative
non_uniq_10 <- lapply(lapply(reads_10,function(i) {i[[1]]==i[[2]]}),unique)
allOne <- sapply(non_uniq_10,length)==1
containTRUE <- sapply(non_uniq_10,'[',1)==T

#length offset of 1 allowed, because filter error might be also 1 base
dupl <- allOne & containTRUE & (reads_length-nachfolger)>=-1 & (reads_length-nachfolger)<=1

fastQ[nrow][dupl==F] -> fastQ

#write output
writeFastq(fastQ,file="output.fastq", mode="w", full=FALSE, compress=F)
