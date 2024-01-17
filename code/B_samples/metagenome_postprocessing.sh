# Generate an index showing which bin each contig is in
for file in $(ls scratch/$1/metabat/*fa); do bin=${file/scratch\/$1\/metabat\/$1./}; bin=${bin/.fa/}; awk -v bin=$bin '/>/{sub(">", "", $0); print $0"\t"bin}' $file; done | sort -V > scratch/$1/bindex.txt

# Generate a list of antismash features by contig
for file in $(ls scratch/$1/antismash/results/NODE*gbk); do short=${file/scratch\/$1\/antismash\/results\//}; short=${short/.region???.gbk/}; printf $short\\t; grep product $file | head -n 1 | grep -Po \".*\" | tr -d "\""; done > scratch/$1/antismash/node_features.txt

# Combine the indices
awk 'FNR==NR{bindex[$1]=$2}; FNR<NR{print $1"\t"bindex[$1]"\t"$2}' scratch/$1/bindex.txt scratch/$1/antismash/node_features.txt | sort -k2,2n -k3 > scratch/$1/antismash/bin_features.txt

# Transform this into a table
cut -f 2-3 scratch/$1/antismash/bin_features.txt | uniq -c | awk '{bins[$2]=1;features[$3]=1;counts[$2,$3]=$1};END{for(bin in bins){printf "\t"bin};print "";for(feature in features){printf feature"\t";for(bin in bins){if((bin,feature) in counts){printf counts[bin,feature]"\t"}else{printf "0\t"}};print ""}}' > scratch/$1/antismash/bin_features_table.txt
