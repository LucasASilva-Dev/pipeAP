# Function: admixture.sh
# License: GPLv3 or later
# Modification date: 22 Mar 2021
# Written by: Marco Milanesi
# Contact: marco.milanesi.mm@gmail.com
# Description: Run Admixture analyses in parallel
# Usage: sh admixture.sh [maxK] [nPar] [ncores] [filename] [email]

#Parameters
maxK=$1 # num of Ks
nPar=$2 #num parallel analyses
ncores=$3 #num of CPU for analyses
filename=$4 # input file
email=$5

#RUN
echo "Start Admixture"
seq 2 $maxK | parallel -j$nPar sh /media/AP_storage/PipeAP/admixture_run.sh ${filename} {} $ncores
echo "END"

#End of the analyses
echo 'Dear,\nthe ADMIXTURE analysis is finished. Check the results!' | mail -s 'Variant calling' $email

