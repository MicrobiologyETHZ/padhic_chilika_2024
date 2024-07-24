#$ -cwd
#$ -S /bin/bash
#$ -N pipeline16s
#$ -V
#$ -pe smp 16
#$ -l h_vmem=2G
#$ -e logs/16s.error.log
#$ -o logs/16s.out.log

# Note that the files have been put into two directories since this was run

output_d="data/processed/16s_both_uparse"

mkdir -p $output_d

rm $output_d/pool.fastq
touch $output_d/pool.fastq

# Label and pool reads
for file in $(ls data/raw/*samples/*.ccs.fastq)
do
    sample=${file/data\/raw\/???_samples\//}
    sample=${sample/.ccs.fastq/}
    sed "s|/ccs.*|/ccs;sample=$sample|" $file
done >> $output_d/pool.fastq

ml cutadapt

cutadapt --discard-untrimmed -g AGRGTTTGATCMTGGCTCAG -a AAGTCGTAACAAGGTAACCC --rc -e 0.1 -m 1000 -o $output_d/pool.primermatch.fasta $output_d/pool.fastq &> $output_d/primermatch.log
sed -i "s| rc$||" $output_d/pool.primermatch.fasta

ml USEARCH

usearch -fastx_uniques $output_d/pool.primermatch.fasta -minuniquesize 1 -sizeout -relabel uniq -fastaout $output_d/uniques.fa -threads 32 &> $output_d/derep.log

usearch -cluster_otus $output_d/uniques.fa -minsize 1 -otus $output_d/otus_uparse.fa -relabel OTU &> $output_d/clustering.log

function lca(){ cat $@ | sed -e '$!{N;s/^\(.*\).*\n\1.*$/\1\n\1/;D;}' | awk -F ";" '{$NF=""; OFS=";"; print $0}'; return; }

usearch -usearch_global $output_d/otus_uparse.fa -db /nfs/cds/Databases/FASTQ_A/SILVA/SILVA138/SILVA_138.1_SSURef_NR99_tax_silva.fasta -id 0.9 -maxaccepts 500 -maxrejects 500 -strand both -top_hits_only -output_no_hits -blast6out $output_d/otus.tax -threads 32 &> $output_d/taxonomy.log

for i in $(cut -f 1 -d $'\t' $output_d/otus.tax | sort | uniq); do id=$(grep -m 1 -P $i'\t' $output_d/otus.tax | cut -f 3 -d$'\t'); res=$(grep -P $i'\t' $output_d/otus.tax | cut -f 2 -d$'\t' | cut -f 1 -d ' ' --complement | lca); echo -e $i'\t'$res'\t'$id; done > $output_d/otus.lca

usearch -otutab $output_d/pool.primermatch.fasta -otus $output_d/otus_uparse.fa -strand both -id 0.97 -otutabout $output_d/otutab.txt -threads 16 &> $output_d/otutab.log
