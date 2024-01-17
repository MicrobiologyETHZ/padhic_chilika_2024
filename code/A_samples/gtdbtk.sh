#$ -cwd
#$ -S /bin/bash
#$ -N GTDBTk
#$ -V
#$ -pe smp 32
#$ -l h_vmem=8G
#$ -e gtdbtk.error.log
#$ -o gtdbtk.out.log

ml GTDBTk

gtdbtk classify_wf --genome_dir gtdbtk_in --out_dir gtdbtk_out --extension fa --cpus 32 --debug
