# Function: seq_variantcalling.R
# License: GPLv3 or later
# Modification date: 11 Mar 2021
# Written by: Marco Milanesi
# Contact: marco.milanesi.mm@gmail.com
# Description: variant calling from BAM file
# Usage: R < seq_alignreads.R  --vanilla --args [refgenome] [bamfile] [regions] [ncores] [outname] [email]

# Libraries ---------------------------------------------------------------------------------------
library(parallel)

# Programs ----------------------------------------------------------------------------------------
bcftools <- "/usr/bin/bcftools"

# User arguments ----------------------------------------------------------------------------------
argv <- commandArgs(T)
refgenome <- argv[1]
bamfile <- argv[2] #for each line a BAM file name with complete PATH
regions <- argv[3] #RefSeq .fai 
ncores <- argv[4]
outname <- argv[5]
email <- argv[6]

# Run Pipeline ------------------------------------------------------------------------------------
#Read regions file
regions <- read.table(file = regions, header = F, stringsAsFactors = F)
regions$NAME <- paste0(regions$V1,":1-",regions$V2)

#Variant calling function
variantcalling <- function(i){
  command <- paste0("bcftools mpileup -C50 -E -a FORMAT/AD,FORMAT/DP,INFO/AD -f ",refgenome," -r ",regions$NAME[i]," -b ",bamfile,
                    " -Ou --threads 1 | bcftools call -Ou -v -m -f GQ,GP --threads 1 | bcftools filter -Oz -g5 -G10 --threads 1 > ",
                    outname,"_",regions$V1[i],".vcf.gz")
  
  system(command)
}

#Run the function
log <- mclapply(X = 1:nrow(regions), FUN = variantcalling, mc.cores = ncores)


#End of analyses - send email
command <- paste0("echo 'Dear,\nthe variant calling is finished. Check the results!' | mail -s 'Variant calling' ",email)
system(command)

