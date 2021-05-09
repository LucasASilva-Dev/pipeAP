# Function: seq_cnvcaller.R
# License: GPLv3 or later
# Modification date: 12 Mar 2021
# Written by: Yuri Tani Utsunomiya
# Contact: ytutsunomiya@gmail.com
# Description: CNVcaller from BAM file 
# Usage: R < seq_cnvcaller.R --vanilla --args [refgenome] [bamfile] [ncores] [email]

# Libraries ---------------------------------------------------------------------------------------
library(parallel)

# Programs ----------------------------------------------------------------------------------------
cnvref <- "/media/AP_storage/AP_storage/yuritani/Projects/NELSEQVAR/CNVcaller/bin/CNVReferenceDB.pl"
cnvrd <- "/media/AP_storage/AP_storage/yuritani/Projects/NELSEQVAR/CNVcaller/Individual.Process.sh"
cnvdiscovery <- "/media/AP_storage/AP_storage/yuritani/Projects/NELSEQVAR/CNVcaller/CNV.Discovery.sh"
cnvgenotype <- "/media/AP_storage/AP_storage/yuritani/Projects/NELSEQVAR/CNVcaller/Genotype.py"
cnvkmer <- "/media/AP_storage/AP_storage/yuritani/Projects/NELSEQVAR/CNVcaller/bin/0.1.Kmer_Generate.py"
cnvlink <- "/media/AP_storage/AP_storage/yuritani/Projects/NELSEQVAR/CNVcaller/bin/0.2.Kmer_Link.py"
perl <- "/usr/bin/perl"
python <- "/usr/bin/python3.6"
blasr <- "/home/yuritani/anaconda3/bin/blasr"
sawriter <- "/home/yuritani/anaconda3/bin/sawriter"

# User arguments ----------------------------------------------------------------------------------
argv <- commandArgs(T)
refgenome <- argv[1]
bamfile <- argv[2]
ncores <- argv[3]
outname <- argv[4]
email <- argv[5]

# Run Pipeline ------------------------------------------------------------------------------------
# Index genome
command <- paste(perl, cnvref, refgenome)
system(command)

# Make genome kmers
command <- paste(python, cnvkmer, refgenome, "800 kmer.fa")
system(command)

#Create .sa file
command <- paste0("mv ",refgenome,".sa ",refgenome,".sa.ORI")
system(command)
command <- paste0(sawriter," ",refgenome,".sa ",refgenome)
system(command)

# Align kmers back to the genome
command <- paste0(blasr," kmer.fa ", refgenome, " --sa ", refgenome, ".sa --out kmer.aln")
command <- paste(command,"-m 5 --noSplitSubreads --minMatch 15 --maxMatch 20 --advanceHalf --advanceExactMatches 10")
command <- paste(command,"--fastMaxInterval --fastSDP --aggressiveIntervalCut --bestn 10 --nproc",ncores)
system(command)

# Make duplicated windows file
command <- paste(python, cnvlink, "kmer.aln 800 window.link")
system(command)

# Parallel RD processing
bamfile <- read.table(file = bamfile, header = F, stringsAsFactors = F)[,1]
rdprocess <- function(file){
  output <- gsub(pattern = "^.+/|\\.bam$", replacement = "", x = file)
  command <- paste("sh",cnvrd,"-b",file,"-h",output,"-d window.link -s X")
  system(command)
}
ncores.par <- min(c(length(bamfile),ncores))
a <- mclapply(X = bamfile, FUN = rdprocess, mc.cores = ncores.par)

# Detect CNVR
rdlist <- list.files(path = "./RD_normalized/")
rdlist <- paste0("./RD_normalized/",rdlist)
write.table(x = rdlist, file = "rdlist.txt", quote = F, row.names = F, col.names = F)
write(x = "", file = "rdexclude.txt")
command <- paste("sh", cnvdiscovery,"-l rdlist.txt -e rdexclude.txt -f 0.1 -h 3 -r 0.5 -p primaryCNVR -m mergeCNVR")
system(command)

# Genotype CNVs
command <- paste(python, cnvgenotype," --cnvfile mergeCNVR --outprefix ",outname," --nproc ",ncores)
system(command)


#End of analyses - send email
command <- paste0("echo 'Dear,\nthe CNVcaller analysis is finished. Check the results!' | mail -s 'Variant calling' ",email)
system(command)

