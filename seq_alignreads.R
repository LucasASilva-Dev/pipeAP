# Function: seq_alignreads.R
# License: GPLv3 or later
# Modification date: 07 Mar 2021
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
dumpfolder <- "/media/TRI_storage/DUMP"

# User arguments ----------------------------------------------------------------------------------
argv <- commandArgs(T)
refgenome <- argv[1]
readfile <- argv[2]
ncores <- argv[3]

# Get readfile ------------------------------------------------------------------------------------
readfile <- read.table(file = readfile, header = T, stringsAsFactors = F)
readfile$out <- gsub(pattern = "\\.bam$", replacement = "", x = readfile$out)

# Run Pipeline ------------------------------------------------------------------------------------
for(i in 1:nrow(readfile)){
  
  # Time stamp
  timestart <- Sys.time()

  # Perform alignments
  command <- paste(bwa, "mem -t", ncores, refgenome, readfile$R1[i], readfile$R2[i])
  command <- paste0(command, " | ", samtools, " view -o ",
                    readfile$out[i], "_unsorted_unmarked.bam -")
  system(command)

  # Sort alignments
  command <- paste0(samtools," sort ", readfile$out[i], "_unsorted_unmarked.bam > ",
                    readfile$out[i], "_unmarked.bam")
  system(command)
  
  # ReadGroup
  name <- basename(readfile$out[i])
  command <- paste0("java -Xmx80g -jar ",picard," AddOrReplaceReadGroups INPUT=",readfile$out[i],"_unmarked.bam",
                    " OUTPUT=",readfile$out[i],"_unmarkedRG.bam"," RGID=",name,
                    " RGLB=",name," RGPL=Illumina RGPU=collapsed RGSM=",name," VALIDATION_STRINGENCY=LENIENT")
  system(command)
  
  # Mark duplicates
  command <- paste0("java -Xmx80g -jar ",picard, " MarkDuplicates ",
                    " VALIDATION_STRINGENCY=LENIENT",
                    " TMP_DIR=", dumpfolder,
                    " ASSUME_SORTED=true REMOVE_DUPLICATES=false",
                    " INPUT=", readfile$out[i], "_unmarkedRG.bam",
                    " METRICS_FILE=", readfile$out[i], "_marked_dup_metrics.txt",
                    " OUTPUT=", readfile$out[i], ".bam")
  system(command)
  
  # Index bam
  command <- paste0(samtools, " index ", readfile$out[i], ".bam")
  system(command)
  
  # Run flagstats
  command <- paste0(samtools, " flagstat ", readfile$out[i], ".bam > ", readfile$out[i], ".flagstat")
  system(command)
  
  # Time stamp
  timeend <- Sys.time()
  timestart
  timeend
  timeend - timestart
  
}
