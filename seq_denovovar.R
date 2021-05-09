# Function: seq_denovovar.R
# License: GPLv3 or later
# Modification date: 24 Mar 2021
# Written by: Marco Milanesi
# Contact: marco.milanesi.mm@gmail.com
# Description: variant calling from BAM file
# Usage: R < seq_denovovar.R  --vanilla --args [vcffile] [trios] [ncores] [email]

# Libraries ---------------------------------------------------------------------------------------
library(parallel)

# Programs ----------------------------------------------------------------------------------------
bcftools <- "/media/AP_storage/AP_storage/SOFTWARE/bcftools-1.12/bin/bcftools-1.12"

# User arguments ----------------------------------------------------------------------------------
argv <- commandArgs(T)
vcffile <- argv[1]
trios <- argv[2]
ncores <- argv[3]
email <- argv[4]

# Run Pipeline ------------------------------------------------------------------------------------

#Function
denovo <- function(x){
  infile <- vcffile[x]
  outfile <- basename(infile)
  outfile <- gsub(pattern = ".vcf.gz", replacement = "", x = outfile)
  
  command <- paste0("bcftools-1.12 +mendelian -m x ",infile," -T ",trios," -Oz -o ",outfile,"_DENOVO.vcf.gz")
  system(command)
  
  command <- paste0("bcftools-1.12 +mendelian ",infile," -T ",trios," -o ",outfile,".log")
  system(command)
}

#Import vcf file
vcffile <- read.table(file = vcffile, header = F, stringsAsFactors = F)[,1]

#Run the function
log <- mclapply(X = 1:length(vcffile), FUN = denovo, mc.cores = ncores)

#End of analyses - send email
command <- paste0("echo 'Dear,\nthe variant calling is finished. Check the results!' | mail -s 'DeNovo variant identification' ",email)
system(command)
