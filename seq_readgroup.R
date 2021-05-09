# Function: seq_readgroup.R
# License: GPLv3 or later
# Modification date: 13 Mar 2021
# Written by: Marco Milanesi
# Contact: marco.milanesi.mm@gmail.com
# Description: picard AddOrReplaceReadGroups
# Usage: R < seq_readgroup.R --vanilla --args [bamfile] [ncores] 

# Libraries ---------------------------------------------------------------------------------------
library(parallel)

# Programs ----------------------------------------------------------------------------------------
picard <- "/usr/bin/picard.jar"

# User arguments ----------------------------------------------------------------------------------
argv <- commandArgs(T)
bamfile <- argv[1]
ncores <- argv[2]

# User functions ----------------------------------------------------------------------------------
rgfun <- function(i){
  outfile <- gsub(pattern = "_old.bam", replacement = ".bam", x = bamfile[i])
  name <- gsub(pattern = "_old.bam", replacement = "", x = basename(bamfile[i]))
  
  command <- paste0("java -jar ",picard," AddOrReplaceReadGroups INPUT=",bamfile[i],
                    " OUTPUT=",outfile," RGID=",name,
                    " RGLB=",name," RGPL=Illumina RGPU=collapsed RGSM=",name," VALIDATION_STRINGENCY=LENIENT")
  system(command)
  
  command <- paste0("samtool index ",outfile)
  system(command)
  
}


# Run Pipeline ------------------------------------------------------------------------------------
bamfile <- read.table(file = bamfile, header = F, stringsAsFactors = F)[,1]

a <- mclapply(X = 1:length(bamfile), FUN = rgfun, mc.cores = ncores)
