# Script to create a file that shows the bin stats for each contig

for sponge in $(ls eggnog)
do
# Cross reference bins with bin stats for each contig
awk 'BEGIN{OFS="\t"};FNR==NR{sub(/_[[:digit:]]+$/,"",$1);var[$1]=$0};FNR!=NR{if(var[$1]!=""){print($2"\t"var[$1])}}' eggnog/$sponge\/$sponge\.sulfur_genes metabat2/$sponge\/prok/$sponge\_bin_assignment.txt > eggnog/$sponge\/$sponge\_bins.sulfur_genes
done
