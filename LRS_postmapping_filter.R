
# #! /bin/Rscript

#### FUNCTIONS

# data reading functions
getDataDelim <- function(pattern, header = T) {
  filenames <- read.table("inSAM_script", header = F, stringsAsFactors = F)[,1]
  names <- sub(pattern, "", filenames)
  names <- sapply(lapply(strsplit(names, "/"), rev), '[', 1)
  tables <- list()
  count <- 0
  for(i in names) {
  	count <- count + 1
    tables[[i]] <- read.delim(filenames[count], stringsAsFactors=F, header= header,comment.char = "", quote = "")
  }
  return(tables)
}

# data writing

writeDataDelim <- function(list) { filenames <- read.table("inSAM_script", header = F, stringsAsFactors = F)[,1]
  								   names <- sub(".sam", "_filt.sam", filenames)
  								   count <- 0
								   for (i in names) {
								   	count <- count + 1
								   	write.table(list[[count]], i, col.names = F, row.names = F, sep = "\t", quote = 1)
								   } }

# filter data
filterData <- function(sam, input) { input_string <- apply(input, 1, function(i) {paste(i[1], i[2], i[3], i[4], sep = "_")})
									 sam_no_match <- lapply(sam, function(sam_table) { sam_string <- apply(sam_table[,1:4], 1, function(i) {paste(i[1], i[2], i[3], i[4], sep = "_")}) 
									 								   no_match <- match(sam_string, input_string, nomatch = 0)==0
									 								   return(sam_table[no_match,]) })
									return(sam_no_match) 
									 								   
									}

#### CODE

input <- read.table("non_unique_assignment", header = F, stringsAsFactors=F)
pattern <- ".sam"
sam <- getDataDelim(pattern, header = F)

sam_no_match <- filterData(sam, input)
writeDataDelim(sam_no_match)
