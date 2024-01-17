#$ -cwd
#$ -S /bin/bash
#$ -N antismash
#$ -V
#$ -pe smp 32
#$ -l h_vmem=2G
#$ -e logs/antismash.error.log
#$ -o logs/antismash.out.log
#$ -t 1-16

i=$(expr $SGE_TASK_ID - 1)
assemblies=($(ls scratch/assemblies/*))
assembly=${assemblies[$i]}
outdir=${assembly/scratch\/assemblies\//}
outdir=${outdir/_spades_scaffolds.fasta/}

echo "Running antismash on $assembly, results in $outdir"

ml antiSMASH

if [ ! -f data/$outdir/index.html ]
then
    antismash -c 32 --output-dir data/$outdir --genefinding-tool prodigal $assembly
fi

