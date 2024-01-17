#$ -cwd
#$ -S /bin/bash
#$ -N emapper
#$ -V
#$ -pe smp 16
#$ -l h_vmem=4G
#$ -e emapper.error.log
#$ -o emapper.out.log
#$ -m a
#$ -M fieldc@ethz.ch
#$ -t 1-5

I=$(expr $SGE_TASK_ID - 1)

DIRS=("20170920.A-C_META1" "20170920.A-MKB_META" "20170920.A-O_META4" "20170920.A-HI_META1" "20170920.A-MKP_META")
DIR=${DIRS[$I]}

mkdir -p eggnog/$DIR

ml eggnog-mapper

emapper.py --dmnd_db /nfs/nas22/fs2201/biol_micro_unix_modules/modules/software/eggnog-mapper/1.0.3-foss-2018b-Python-2.7.15/lib/python2.7/site-packages/eggnogmapper-1.0-py2.7.egg/data/eggnog_proteins.dmnd -o eggnog/$DIR/$DIR -i annotation/kaiju/$DIR/$DIR\.faa --cpu 16 -m diamond
