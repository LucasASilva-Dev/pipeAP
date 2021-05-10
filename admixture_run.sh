filename=$1
K=$2
ncpu=$3

admixture --cv ${filename}".bed" $K -j$ncpu | tee ${filename}.${K}.log

