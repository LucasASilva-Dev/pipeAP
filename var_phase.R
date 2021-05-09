# Function: var_phase.R
# License: GPLv3 or later
# Modification date: 12 Mar 2021
# Written by: Yuri Tani Utsunomiya
# Contact: ytutsunomiya@gmail.com
# Description: phase variant data
# Usage: Rscript var_phase.R [bedfile] [autonum] [outfile] [ncores]

# Programs ----------------------------------------------------------------------------------------
plink <- "/usr/bin/plink_v1.90b5.4"
eagle <- "/usr/bin/eagle"
tabix <- "/usr/bin/tabix"

# User arguments ----------------------------------------------------------------------------------
argv <- commandArgs(T)
bedfile <- argv[1]
autonum <- as.integer(argv[2])
outfile <- argv[3]
ncores <- argv[4]

# Temporary file ID generator ---------------------------------------------------------------------
tmpfile <- paste0("tmp_", paste(sample(LETTERS, size = 10, replace = T), collapse=""))
tmpfile <- gsub(pattern = "\\.", replacement = "", paste0(tmpfile,as.numeric(Sys.time())))

# Build genetic map ---------------------------------------------------------------------------------------------
bim <- read.table(file = paste0(bedfile,".bim"), header = F, stringsAsFactors = F)
genmap <- bim[,c(1,4,4,4)]
genmap[,3] <- 1
genmap[,4] <- genmap[,4]/1e+6
colnames(genmap) <- c("chr","position","COMBINED_rate(cM/Mb)","Genetic_Map(cM)")

# Run phasing algorithm ---------------------------------------------------------------------------
for(i in 1:autonum){
  
  # Make VCF file
  command <- paste(plink, "--chr-set", autonum,
                   "--bfile", bedfile,
                   "--chr",i,
                   "--recode vcf-iid",
                   "--out", tmpfile)
  system(command)
  
  # Index with tabix
  command <- paste0("cat ", tmpfile, ".vcf | bgzip > ", tmpfile, ".vcf.gz")
  system(command)
  command <- paste0(tabix, " -p vcf ", tmpfile, ".vcf.gz")
  system(command)
  
  # Export genetic map
  write.table(x = genmap[which(genmap$chr == i),], file = paste0(tmpfile, "_genmap.txt"), quote = F,
              sep = " ", row.names = F, col.names = T)
  
  # Phasing
  command <- paste0(eagle, " --vcf ", tmpfile, ".vcf.gz",
                   " --geneticMapFile=", tmpfile, "_genmap.txt",
                   " --chromX ", autonum+1,
                   " --numThreads ", ncores,
                   " --outPrefix=", outfile, "_chr", i)
  system(command)
  
  # Clear directory
  command <- paste0("rm ", tmpfile, "*")
  system(command)
  
}