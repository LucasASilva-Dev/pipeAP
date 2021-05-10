# Function: seq_alignreads.R
# License: GPLv3 or later
# Modification date: 04 Mar 2021
# Written by: Yuri Tani Utsunomiya
# Contact: ytutsunomiya@gmail.com
# Description: read alignment against reference genome
# Usage: Rscript seq_alignreads.R [refgenome] [readsfile] [ncores]

# Libraries ---------------------------------------------------------------------------------------
library(parallel)

# Programs ----------------------------------------------------------------------------------------
bwa <- "/usr/bin/bwa"
samtools <- "/usr/bin/samtools"
picard <- "/usr/bin/picard.jar"

# User arguments ----------------------------------------------------------------------------------
argv <- commandArgs(T)
refgenome <- argv[1]
readfile <- argv[2]
ncores <- argv[3]

# Time stamp --------------------------------------------------------------------------------------
timestart <- Sys.time()

# Get readfile ------------------------------------------------------------------------------------
readfile <- read.table(file = readfile, header = T, stringsAsFactors = F)
readfile$out <- gsub(pattern = "\\.bam$", replacement = "", x = readfile$out)

# Perform alignments ------------------------------------------------------------------------------
for(i in 1:nrow(readfile)){
  command <- paste(bwa, "mem -t", ncores, refgenome, readfile$R1[i], readfile$R2[i])
  command <- paste0(command, " | ", samtools, " view -o ",
                    readfile$out[i], "_unsorted_unmarked.bam -")
  system(command)
}

# Sort alignments ----------------------------------------------------------------------------------
sortalign <- function(i){
  command <- paste0(samtools," sort ", readfile$out[i], "_unsorted_unmarked.bam > ",
                    readfile$out[i], "_unmarked.bam")
  system(command)
}
log <- mclapply(X = 1:nrow(readfile), FUN = sortalign, mc.cores = ncores)

# Mark duplicates ----------------------------------------------------------------------------------
markdup <- function(i){
  command <- paste0("java -Xmx10g -jar ",picard, " MarkDuplicates ",
                    " VALIDATION_STRINGENCY=LENIENT",
                    " ASSUME_SORTED=true REMOVE_DUPLICATES=false",
                    " INPUT=", readfile$out[i], "_unmarked.bam",
                    " METRICS_FILE=", readfile$out[i], "_marked_dup_metrics.txt",
                    " OUTPUT=", readfile$out[i], ".bam")
  system(command)
}
log <- mclapply(X = 1:nrow(readfile), FUN = markdup, mc.cores = ncores)

# Index bam ----------------------------------------------------------------------------------------
indexbam <- function(i){
  command <- paste0(samtools, " index ", readfile$out[i], ".bam")
  system(command)
}
log <- mclapply(X = 1:nrow(readfile), FUN = indexbam, mc.cores = ncores)

# Run flagstats ------------------------------------------------------------------------------------
flagstats <- function(i){
  command <- paste0(samtools, " flagstat ", readfile$out[i], ".bam > ", readfile$out[i], ".flagstat")
  system(command)
}
log <- mclapply(X = 1:nrow(readfile), FUN = flagstats, mc.cores = ncores)

# Remove temporary files ---------------------------------------------------------------------------
for(i in 1:nrow(readfile)){
  command <- paste0("rm ", readfile$out[i], "_*unmarked.bam")
  system(command)
}

# Time stamp --------------------------------------------------------------------------------------
timeend <- Sys.time()
timestart
timeend
timeend - timestart
