# Script to create a file that shows the bin stats for each contig

for sponge in $(ls raw)
do

# Assign bin number to each contig
for file in $(ls metabat2/$sponge\/prok/bins/$sponge\_prok.*); do n=${file/metabat2\/$sponge\/prok\/bins\/$sponge\_prok./}; n=${n/.fa/}; awk -v n=${n} -v s=${sponge} '/>/{sub(">","",$0);print $0"\t"s"_prok."n}' $file;done | sort -V > metabat2/$sponge\/prok/$sponge\_bin_assignment.txt

# Cross reference bins with bin stats for each contig
awk 'FNR==NR{var[$1]=$1"\t"$13"\t"$14"\t"$15} FNR!=NR{print($1"\t"var[$2])}' checkm/$sponge\/prok/$sponge\_prok_checkm.stats metabat2/$sponge\/prok/$sponge\_bin_assignment.txt > checkm/$sponge\/prok/$sponge\_prok_contig_checkm.stats

done
