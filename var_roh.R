# Function: var_roh.R
# License: GPLv3 or later
# Modification date: 12 Mar 2021
# Written by: Yuri Tani Utsunomiya
# Contact: ytutsunomiya@gmail.com
# Description: detect runs of homozygosity from variant data
# Usage: Rscript var_roh.R [bedfile] [autonum] [outfile] [ncores]

# Libraries ---------------------------------------------------------------------------------------
library(RZooRoH)
library(parallel)
library(data.table)

# Programs ----------------------------------------------------------------------------------------
plink <- "/usr/bin/plink_v1.90b5.4"

# User arguments ----------------------------------------------------------------------------------
argv <- commandArgs(T)
bedfile <- argv[1]
autonum <- argv[2]
outfile <- argv[3]
ncores <- argv[4]

# Temporary file ID generator ---------------------------------------------------------------------
tmpfile <- paste0("tmp_", paste(sample(LETTERS, size = 10, replace = T), collapse=""))
tmpfile <- gsub(pattern = "\\.", replacement = "", paste0(tmpfile,as.numeric(Sys.time())))

# Format input ------------------------------------------------------------------------------------
command <- paste(plink, "--autosome --chr-set", autonum,
                 "--bfile", bedfile,
                 "--recode oxford",
                 "--out", tmpfile)
system(command)

# Get genome size ---------------------------------------------------------------------------------
command <- paste0("cut -d' ' -f1,3 ", tmpfile, ".gen > ", tmpfile,"_genomesize.txt")
system(command)
genomesize <- fread(paste0(tmpfile,"_genomesize.txt"), head=F, sep=" ", stringsAsFactors = F)
mrkdist <- diff(genomesize$V2)
mrkdist <- mrkdist[which(mrkdist > 0)]
genomesize <- sum(as.numeric(mrkdist))

# Load gen file -----------------------------------------------------------------------------------
command <- paste0("tail -n+3 ", tmpfile, ".sample | cut -d' ' -f2 > ", tmpfile, "_ids.txt")
system(command)
geno <- zoodata(paste0(tmpfile, ".gen"), zformat = "gp", chrcol=1,
                poscol=3, supcol=5, samplefile = paste0(tmpfile, "_ids.txt"))

# Detect ROH --------------------------------------------------------------------------------------
model <- zoomodel(K = 2)
results <- zoorun(zoomodel = model, zooin = geno, nT = ncores)
results@hbdseg$id <- geno@sample_ids[results@hbdseg$id]

# Compute FROH ------------------------------------------------------------------------------------
getfroh <- function(result, genomesize){
  auto <- as.data.frame(matrix(nrow = length(unique(result@hbdseg$id)),ncol = 16,NA))
  names(auto) <- c("ID",
                   "FROH1", "FROH2", "FROH4", "FROH8", "FROH16",
                   "MROH1", "MROH2", "MROH4", "MROH8", "MROH16",
                   "NROH1", "NROH2", "NROH4", "NROH8", "NROH16")
  auto$ID <- sort(unique(result@hbdseg$id))
  rohsize <- c(1, 2, 4, 8, 16)*1000000
  for (i in 1:nrow(auto)){
    a <- result@hbdseg[which(result@hbdseg$id==auto$ID[i]),]
    for (j in 1:length(rohsize)){
      segs <- which(a$length>rohsize[j])
      auto[i,paste0("FROH",rohsize[j]/1000000)] <- sum(a[segs,"length"])/genomesize
      auto[i,paste0("MROH",rohsize[j]/1000000)] <- median(a[segs,"length"])
      auto[i,paste0("NROH",rohsize[j]/1000000)] <- length(segs)
    }
  }
  return(auto)
}
froh <- getfroh(results, genomesize)

# Make output ------------------------------------------------------------------------------------
roh <- results@hbdseg[,c("id","chrom","start_pos","end_pos","length")]
colnames(roh) <- c("ID","CHR","BP1","BP2","LENGTH")
write.table(x = roh, file = paste0(outfile,".roh"),
            append = F, quote = F, row.names = F, col.names = T)
write.table(x = froh, file = paste0(outfile,".froh"),
            append = F, quote = F, row.names = F, col.names = T)

# Clear temporary files --------------------------------------------------------------------------
command <- paste0("rm ", tmpfile, "*")
system(command)

