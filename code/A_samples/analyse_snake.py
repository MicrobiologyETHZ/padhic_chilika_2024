import glob

ROOTFOLDER = '/nfs/cds-shini.ethz.ch/exports/biol_micro_cds_gr_sunagawa/scratch/hans/sponges'
RAWFOLDER = ROOTFOLDER + '/raw/'
RENAMEDFOLDER = ROOTFOLDER + '/renamed/'
QCFOLDER = ROOTFOLDER + '/qc/'
ASSEMBLYFOLDER = ROOTFOLDER + '/assembly/'
METABATFOLDER = ROOTFOLDER + '/metabat2/'
CHECKMFOLDER = ROOTFOLDER + '/checkm/'
EUKREPFOLDER = ROOTFOLDER + '/eukrep/'
ANNOFOLDER = ROOTFOLDER + '/annotation/'

RAW_FILES = glob.glob(RAWFOLDER + '*/*.fastq.gz')
RENAMED_FILES = []
PAIRING_FILES = set()
QC_FILES = set()
ASSEMBLYFILES = []
FASTQC_FILES = set()
COUNT_FILES = []
ASSEMBLY_FILES = set()
LONG_CONTIG_FILES = set()
CONTIG_BWA_INDEX_FILES = set()
CONTIG_ALIGNMENT_FILES = set()
METABAT_MARKER_FILES = set()
CHECKM_MARKER_FILES = set()
EUKREP_FILES = set()
ANNO_FILES = set()


for fqf in RAW_FILES:
    renamed = fqf.replace(RAWFOLDER, RENAMEDFOLDER).replace('_R','.').replace('fastq','fq')
    pairing = fqf.replace(RAWFOLDER, RENAMEDFOLDER).split('_R')[0] + '.pairing.txt'
    qc = fqf.replace(RAWFOLDER, QCFOLDER).replace('_R','.').replace('fastq','fq')
    qc_single = fqf.replace(RAWFOLDER, QCFOLDER).split('_R')[0] + '.SINGLETONS.fq.gz'
    fqc_renamed = renamed.replace('fq.gz', 'fastqc.done')
    fqc_qc = qc.replace('fq.gz', 'fastqc.done')
    fqc_qc_single = qc_single.replace('fq.gz', 'fastqc.done')
    assembly = fqf.replace(RAWFOLDER, ASSEMBLYFOLDER).split('_R')[0] + '.spades.done'
    longcontigs = fqf.replace(RAWFOLDER, ASSEMBLYFOLDER).split('_R')[0] + '_spades/contigs_min5K.fasta'
    
    
    bwaalign_prok = os.path.dirname(fqf.replace(RAWFOLDER, METABATFOLDER)) + '/prok/' + os.path.basename(fqf.replace(RAWFOLDER, METABATFOLDER)).split('_R')[0]  + '.bam'
    bwaalign_prok_marker = bwaalign_prok.replace('.bam', '.done')
    
    bwaalign_euk = os.path.dirname(fqf.replace(RAWFOLDER, METABATFOLDER)) + '/euk/' + os.path.basename(fqf.replace(RAWFOLDER, METABATFOLDER)).split('_R')[0]  + '.bam'
    bwaalign_euk_marker = bwaalign_euk.replace('.bam', '.done')
    
    bwaalign_all = os.path.dirname(fqf.replace(RAWFOLDER, METABATFOLDER)) + '/all/' + os.path.basename(fqf.replace(RAWFOLDER, METABATFOLDER)).split('_R')[0]  + '.bam'
    bwaalign_all_marker = bwaalign_all.replace('.bam', '.done')

    metabat_prok = os.path.dirname(fqf.replace(RAWFOLDER, METABATFOLDER)) + '/prok/' + os.path.basename(fqf.replace(RAWFOLDER, METABATFOLDER)).split('_R')[0]  + '_metabat.done'
    metabat_euk = os.path.dirname(fqf.replace(RAWFOLDER, METABATFOLDER)) + '/euk/' + os.path.basename(fqf.replace(RAWFOLDER, METABATFOLDER)).split('_R')[0]  + '_metabat.done'
    metabat_all = os.path.dirname(fqf.replace(RAWFOLDER, METABATFOLDER)) + '/all/' + os.path.basename(fqf.replace(RAWFOLDER, METABATFOLDER)).split('_R')[0]  + '_metabat.done' 
    
    checkm_euk = os.path.dirname(fqf.replace(RAWFOLDER, CHECKMFOLDER)) + '/euk/' + os.path.basename(fqf.replace(RAWFOLDER, CHECKMFOLDER)).split('_R')[0]  + '_checkm.done'
    checkm_prok = os.path.dirname(fqf.replace(RAWFOLDER, CHECKMFOLDER)) + '/prok/' + os.path.basename(fqf.replace(RAWFOLDER, CHECKMFOLDER)).split('_R')[0]  + '_checkm.done'    
    checkm_all = os.path.dirname(fqf.replace(RAWFOLDER, CHECKMFOLDER)) + '/all/' + os.path.basename(fqf.replace(RAWFOLDER, CHECKMFOLDER)).split('_R')[0]  + '_checkm.done'

    checkm_euk_ssu = os.path.dirname(fqf.replace(RAWFOLDER, CHECKMFOLDER)) + '/euk/' + os.path.basename(fqf.replace(RAWFOLDER, CHECKMFOLDER)).split('_R')[0]  + '_checkm_ssu_finder.done'
    checkm_prok_ssu = os.path.dirname(fqf.replace(RAWFOLDER, CHECKMFOLDER)) + '/prok/' + os.path.basename(fqf.replace(RAWFOLDER, CHECKMFOLDER)).split('_R')[0]  + '_checkm_ssu_finder.done'
    checkm_all_ssu = checkm_prok = os.path.dirname(fqf.replace(RAWFOLDER, CHECKMFOLDER)) + '/all/' + os.path.basename(fqf.replace(RAWFOLDER, CHECKMFOLDER)).split('_R')[0]  + '_checkm_ssu_finder.done'
    euk = os.path.dirname(fqf.replace(RAWFOLDER, EUKREPFOLDER)) + '/' + os.path.basename(fqf.replace(RAWFOLDER, EUKREPFOLDER)).split('_R')[0]  + '_euk_contigs_min5K.fasta'
    prok = os.path.dirname(fqf.replace(RAWFOLDER, EUKREPFOLDER)) + '/' + os.path.basename(fqf.replace(RAWFOLDER, EUKREPFOLDER)).split('_R')[0]  + '_prok_contigs_min5K.fasta'
   
    prodigal_all_marker = os.path.dirname(fqf.replace(RAWFOLDER, ANNOFOLDER)) + '/' + os.path.basename(fqf.replace(RAWFOLDER, ANNOFOLDER)).split('_R')[0]  + '_prodigal.done'
    prodigal_all_fna = os.path.dirname(fqf.replace(RAWFOLDER, ANNOFOLDER)) + '/' + os.path.basename(fqf.replace(RAWFOLDER, ANNOFOLDER)).split('_R')[0]  + '.fna'
    prodigal_all_faa = os.path.dirname(fqf.replace(RAWFOLDER, ANNOFOLDER)) + '/' + os.path.basename(fqf.replace(RAWFOLDER, ANNOFOLDER)).split('_R')[0]  + '.faa'
    prodigal_all_gff = os.path.dirname(fqf.replace(RAWFOLDER, ANNOFOLDER)) + '/' + os.path.basename(fqf.replace(RAWFOLDER, ANNOFOLDER)).split('_R')[0]  + '.gff'
    
    allcontigs = os.path.dirname(fqf.replace(RAWFOLDER, EUKREPFOLDER)) + '/' + os.path.basename(fqf.replace(RAWFOLDER, EUKREPFOLDER)).split('_R')[0]  + '_all_contigs_min5K.fasta'
    
    RENAMED_FILES.append(renamed)
    QC_FILES.add(qc)
    QC_FILES.add(qc_single)
    FASTQC_FILES.add(fqc_renamed)
    FASTQC_FILES.add(fqc_qc)
    FASTQC_FILES.add(fqc_qc_single)
    PAIRING_FILES.add(pairing)
    ASSEMBLY_FILES.add(assembly)
    LONG_CONTIG_FILES.add(longcontigs)
    CONTIG_ALIGNMENT_FILES.add(bwaalign_prok)
    CONTIG_ALIGNMENT_FILES.add(bwaalign_prok_marker)
    CONTIG_ALIGNMENT_FILES.add(bwaalign_euk)
    CONTIG_ALIGNMENT_FILES.add(bwaalign_euk_marker)
    CONTIG_ALIGNMENT_FILES.add(bwaalign_all)
    CONTIG_ALIGNMENT_FILES.add(bwaalign_all_marker)
    METABAT_MARKER_FILES.add(metabat_prok)
    METABAT_MARKER_FILES.add(metabat_euk)
    METABAT_MARKER_FILES.add(metabat_all)
    CHECKM_MARKER_FILES.add(checkm_euk)
    CHECKM_MARKER_FILES.add(checkm_prok)
    CHECKM_MARKER_FILES.add(checkm_all)
    CHECKM_MARKER_FILES.add(checkm_euk_ssu)
    CHECKM_MARKER_FILES.add(checkm_all_ssu)
    CHECKM_MARKER_FILES.add(checkm_prok_ssu)


    EUKREP_FILES.add(euk)
    EUKREP_FILES.add(prok)
    EUKREP_FILES.add(allcontigs)
    
    ANNO_FILES.add(prodigal_all_marker)
    ANNO_FILES.add(prodigal_all_fna)
    ANNO_FILES.add(prodigal_all_faa)
    ANNO_FILES.add(prodigal_all_gff)
QC_FILES = list(QC_FILES)
FASTQC_FILES = list(FASTQC_FILES)
PAIRING_FILES = list(PAIRING_FILES)
ASSEMBLY_FILES = list(ASSEMBLY_FILES)
LONG_CONTIG_FILES = list(LONG_CONTIG_FILES)
CONTIG_BWA_INDEX_FILES = list(CONTIG_BWA_INDEX_FILES)
CONTIG_ALIGNMENT_FILES = list(CONTIG_ALIGNMENT_FILES)
METABAT_MARKER_FILES = list(METABAT_MARKER_FILES)
CHECKM_MARKER_FILES = list(CHECKM_MARKER_FILES)
EUKREP_FILES = list(EUKREP_FILES)
ANNO_FILES = list(ANNO_FILES)
rule all:
    input:
        RENAMED_FILES,
        PAIRING_FILES,
        QC_FILES,
        ASSEMBLY_FILES,
        LONG_CONTIG_FILES,
        CONTIG_ALIGNMENT_FILES,
        METABAT_MARKER_FILES,
        CHECKM_MARKER_FILES,
        EUKREP_FILES
        ANNO_FILES
        
rule rename:
    input:
        fq = RAWFOLDER + '{SAMPLE}_R{PAIR}.fastq.gz'
    output:
        fq = RENAMEDFOLDER + '{SAMPLE}.{PAIR}.fq.gz'
    params:
        lsferrfile = RENAMEDFOLDER + '{SAMPLE}.{PAIR}.rename.lsferr',
        lsfoutfile = RENAMEDFOLDER + '{SAMPLE}.{PAIR}.rename.lsfout',
        scratch = 10000,
        mem = 2000,
        time = 100
    benchmark:
        RENAMEDFOLDER + '{SAMPLE}.{PAIR}.rename.benchmark'
    threads:
        1
    shell:
        'gunzip -c {input.fq} | python rename.py | gzip -1 -c > {output.fq}'


rule reformat_vpair:
    input:
        fq1 = '{sample}.1.fq.gz',
        fq2 = '{sample}.2.fq.gz',
    output:
        stats = '{sample}.pairing.txt'
    params:
        lsfoutfile = '{sample}.pairing.lsfout.log',
        lsferrfile = '{sample}.pairing.lsferr.log',
        scratch = 1000,
        mem = 4000,
        time = 235
    threads:
        1
    shell:
        'reformat.sh -Xmx1G threads={threads} qin=33 in1={input.fq1} in2={input.fq2} vpair 2> {output.stats} || echo $?'

rule qc:
    input:
        fqgz1 = RENAMEDFOLDER + '{sample}.1.fq.gz',
        fqgz2 = RENAMEDFOLDER + '{sample}.2.fq.gz'
    output:
        fqgz1 = QCFOLDER + '{sample}.1.fq.gz',
        fqgz2 = QCFOLDER + '{sample}.2.fq.gz',
        fqgz_qc_failed = QCFOLDER + '{sample}_qc_failed.fastq.gz',
        fqgz_qc_singletons = QCFOLDER + '{sample}.SINGLETONS.fq.gz',
        qc_stats = QCFOLDER + '{sample}_qc.stats',
        fqgz_adapter_matched = QCFOLDER + '{sample}_adapter.matched.fastq.gz',
        fqgz_adapter_singletons = QCFOLDER + '{sample}_adapter.singletons.fastq.gz',
        adapter_stats = QCFOLDER + '{sample}_adapter.stats',
        fqgz_phix_matched = QCFOLDER + '{sample}_phix.matched.fastq.gz',
        fqgz_phix_singletons = QCFOLDER + '{sample}_phix.singletons.fastq.gz',
        phix_stats = QCFOLDER + '{sample}_phix.stats',
    params:
        lsferrfile = QCFOLDER + '{sample}_qc_lsferr.log',
        lsfoutfile = QCFOLDER + '{sample}_qc_lsfout.log',
        scratch = 1000,
        mem = 1000,
        time = 100,
    log:
        qc_log = QCFOLDER + '{sample}_qc.log',
        adapter_log = QCFOLDER + '{sample}_adapter.log',
        phix_log = QCFOLDER + '{sample}_phix.log'
    threads:
        32
    benchmark:
        QCFOLDER + '{sample}.qc.benchmark'
    shell:
        'bbduk.sh -Xmx1G usejni=t threads=14 overwrite=t qin=33 in1={input.fqgz1} in2={input.fqgz2} ref=/nfs/cds-shini.ethz.ch/exports/biol_micro_cds_gr_sunagawa/SequenceStorage/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 out=stdout.fq outm={output.fqgz_adapter_matched} outs={output.fqgz_adapter_singletons} refstats={output.adapter_stats} statscolumns=5 2> {log.adapter_log} | bbduk.sh -Xmx1G usejni=t threads=4 interleaved=true overwrite=t qin=33 in=stdin.fq out=stdout.fq outm={output.fqgz_phix_matched} outs={output.fqgz_phix_singletons} ref=/nfs/cds-shini.ethz.ch/exports/biol_micro_cds_gr_sunagawa/SequenceStorage/resources/phix174_ill.ref.fa.gz k=31 hdist=1 refstats={output.phix_stats} statscolumns=5 2> {log.phix_log} | bbduk.sh -Xmx1G usejni=t threads=14 overwrite=t interleaved=true qin=33 in=stdin.fq fastawrap=10000 out1={output.fqgz1} out2={output.fqgz2} interleaved=true outm={output.fqgz_qc_failed} outs={output.fqgz_qc_singletons} minlength=45 maq=20 maxns=0 overwrite=t stats={output.qc_stats} statscolumns=5 trimq=14 qtrim=r 2> {log.qc_log}'


import os
def getMetaspdesOutputFolder(wildcards): 
    return ASSEMBLYFOLDER + wildcards.sample + '_spades'
def getMetaspdesOutputFolderName(wildcards):
    return os.path.basename(getMetaspdesOutputFolder(wildcards))


def getMetaSpadesParams(sample):
    infiles = [QCFOLDER + sample + '.1.fq.gz', QCFOLDER + sample +'.2.fq.gz', QCFOLDER + sample +'.SINGLETONS.fq.gz']
    totalSize = 0 #in megabytes
    for infile in infiles:
        totalSize = totalSize + (os.path.getsize(infile) >> 20)
    minMem = 64000
    minScratch = 64000
    mem = totalSize * 2 * 7
    scratch = totalSize * 2 * 7
    if mem < minMem:
        mem = minMem
    if scratch < minScratch:
        scratch = minScratch
    spadesMem = int((mem*0.95)/1000)
    data =  {'time': '7150', 'mem': str(mem), 'spadesMem': str(spadesMem), 'scratch': str(scratch) }
    return data

def getFq1Name(wildcards):
    return os.path.basename(QCFOLDER + wildcards.sample + '.1.fq.gz')
def getFq2Name(wildcards):
    return os.path.basename(QCFOLDER + wildcards.sample + '.2.fq.gz')
def getFqsName(wildcards):
    return os.path.basename(QCFOLDER + wildcards.sample + '.SINGLETONS.fq.gz')
def getMetaSpadesMem(wildcards):
    return getMetaSpadesParams(wildcards.sample)['mem']
def getMetaSpadesScratch(wildcards):
    return getMetaSpadesParams(wildcards.sample)['scratch']
def getMetaSpadesSpadesMem(wildcards):
    return getMetaSpadesParams(wildcards.sample)['spadesMem']
def getMetaSpadesTime(wildcards):
    return getMetaSpadesParams(wildcards.sample)['time']
def test(wildcards):
    s = wildcards.sample
    print(s)
    return ''

rule assemble_metaspades:
    input:
        fq1 = QCFOLDER + '{sample}.1.fq.gz',
        fq2 = QCFOLDER + '{sample}.2.fq.gz',
        fqs = QCFOLDER + '{sample}.SINGLETONS.fq.gz'
    output:
        marker = touch(ASSEMBLYFOLDER + '{sample}.spades.done')
    params: 
        lsferrfile = ASSEMBLYFOLDER + '{sample}.spades.lsferr.log',
        lsfoutfile = ASSEMBLYFOLDER + '{sample}.spades.lsfout.log',
        fq1 = getFq1Name,
        fq2 = getFq2Name,
        fqs = getFqsName,
        time = 7150, #getMetaSpadesTime,
        spadesMem = 200,#getMetaSpadesSpadesMem,
        mem = 220000,#getMetaSpadesMem,
        scratch = 100000,#getMetaSpadesScratch,
        outputfolder = getMetaspdesOutputFolder,
        outputfoldername = getMetaspdesOutputFolderName,
    benchmark:
        ASSEMBLYFOLDER + '{sample}.spades.benchmark'
    log:
        ASSEMBLYFOLDER + '{sample}.spades.log'
    threads:
        32
    shell:
        'metaspades.py -t {threads} -m 2900 --only-assembler --pe1-1 {input.fq1} --pe1-2 {input.fq2} --pe1-s {input.fqs} -o {params.outputfolder} &> {log}'


rule remove_short_contigs:
    input:
        spades_done = ASSEMBLYFOLDER + '{sample}/{sample}.spades.done'
    output:
        contigs = ASSEMBLYFOLDER + '{sample}/{sample}_spades/contigs_min5K.fasta'
    params:
        contigs = ASSEMBLYFOLDER + '{sample}/{sample}_spades/contigs.fasta'
    benchmark:
        ASSEMBLYFOLDER + '{sample}.remove_short_contigs.benchmark'
    shell:
        '''
        /nfs/home/hansr/miniconda3/envs/biopython/bin/python /nfs/home/hansr/shortSequenceRemover.py {params.contigs} 5000 > {output.contigs}
        '''


        '''
        cp {input.contigs} {params.outputpath}
        bwa index {params.outputpath} 2> {log}
        rm {params.outputpath}
        '''

rule bwa_align:
    input:
        fq1 = QCFOLDER + '{sample}/{sample}.1.fq.gz',
        fq2 = QCFOLDER + '{sample}/{sample}.2.fq.gz',
        fqs = QCFOLDER + '{sample}/{sample}.SINGLETONS.fq.gz',
        contigs = EUKREPFOLDER + '{sample}/{sample}_{euk_prok}_contigs_min5K.fasta'
    output:
        bam = METABATFOLDER + '{sample}/{euk_prok}/{sample}.bam',
        marker = touch(METABATFOLDER + '{sample}/{euk_prok}/{sample}.done')
    params:
        bwapath = METABATFOLDER + '{sample}/{euk_prok}/{sample}_{euk_prok}_contigs_min5K.fasta',
        sam = METABATFOLDER + '{sample}/{euk_prok}/{sample}.sam'
    log:
        METABATFOLDER + '{sample}/{euk_prok}/bwa.align.log'
    benchmark:
        METABATFOLDER + '{sample}/{euk_prok}/bwa.align.benchmark'
    threads:
        32
    shell:
        '''
        cp {input.contigs} {params.bwapath}
        bwa index {params.bwapath} &>> {log}
        bwa mem -t {threads} {params.bwapath} {input.fq1} {input.fq2} 2>> {log}| samtools view -h -F 4 - | python /nfs/cds/mOTUv2/mOTUs_v2/bin/msamtools_python.py 0.97 200 90 None > {params.sam}
        bwa mem -t {threads} {params.bwapath} {input.fqs} 2>> {log}| samtools view -F 4 - | python /nfs/cds/mOTUv2/mOTUs_v2/bin/msamtools_python.py 0.97 200 90 None >> {params.sam}
        samtools view -bu {params.sam} | samtools sort -@ {threads} - > {output.bam} 2>> {log}
        samtools index {output.bam} 2>> {log}
        
        '''

rule metabat:
    input:
        bam = METABATFOLDER + '{sample}/{euk_prok}/{sample}.bam',
        contigs = EUKREPFOLDER + '{sample}/{sample}_{euk_prok}_contigs_min5K.fasta',
        marker = METABATFOLDER + '{sample}/{euk_prok}/{sample}.done'
    output:
        marker = touch(METABATFOLDER + '{sample}/{euk_prok}/{sample}_metabat.done')
    params:
        metabatbinfolder = METABATFOLDER + '{sample}/{euk_prok}/bins/{sample}_{euk_prok}',
        depth = METABATFOLDER + '{sample}/{euk_prok}/{sample}_depth.txt'
    log:
        METABATFOLDER + '{sample}/{euk_prok}/metabat.log'
    benchmark:
        METABATFOLDER + '{sample}/{euk_prok}/metabat.benchmark'
    threads:
        16
    shell:
        '''
        /nfs/practical/BlockCourse_Fall2017/week2/resources/metabat/jgi_summarize_bam_contig_depths --outputDepth {params.depth} {input.bam} &>> {log}
        /nfs/practical/BlockCourse_Fall2017/week2/resources/metabat/metabat2 -t {threads} -i {input.contigs} -o {params.metabatbinfolder} --unbinned -v -a {params.depth} &>> {log}
        '''

rule checkm:
    input:
        marker = METABATFOLDER + '{sample}/{euk_prok}/{sample}_metabat.done',
        bam = METABATFOLDER + '{sample}/{euk_prok}/{sample}.bam',
        contigs = EUKREPFOLDER + '{sample}/{sample}_{euk_prok}_contigs_min5K.fasta'
    output:
        marker = touch(CHECKMFOLDER + '{sample}/{euk_prok}/{sample}_checkm.done')
    params:
        tetra = CHECKMFOLDER + '{sample}/{euk_prok}/{sample}_{euk_prok}_tetra.tsv',
        coverage = CHECKMFOLDER + '{sample}/{euk_prok}/{sample}_{euk_prok}_coverage.tsv',
        coverage_profile = CHECKMFOLDER + '{sample}/{euk_prok}/{sample}_{euk_prok}_coverage_profile.txt',
        bins = METABATFOLDER + '{sample}/{euk_prok}/bins/',
        chmfolder = CHECKMFOLDER + '{sample}/{euk_prok}/checkm/',
        plotsfolder = CHECKMFOLDER + '{sample}/{euk_prok}/plots/',
        stats = CHECKMFOLDER + '{sample}/{euk_prok}/{sample}_{euk_prok}_checkm.stats'
    log:
        CHECKMFOLDER + '{sample}/{euk_prok}/checkm.log'
    benchmark:
        CHECKMFOLDER + '{sample}/{euk_prok}/checkm.benchmark'
    threads:
        3
    shell:
        '''
        checkm lineage_wf -f {params.stats} -t {threads} -x fa --pplacer_threads {threads} {params.bins} {params.chmfolder} 2> {log}
        checkm bin_qa_plot -x fa --image_type pdf {params.chmfolder} {params.bins} {params.plotsfolder} 2>> {log}
        #checkm bin_qa_plot -x fa --image_type png {params.chmfolder} {params.bins} {params.plotsfolder} 2>> {log}
        checkm tetra {input.contigs} {params.tetra} 2>> {log}
        checkm dist_plot -x fa --image_type pdf {params.chmfolder} {params.bins} {params.plotsfolder} {params.tetra} 95 2>> {log}
        #checkm dist_plot -x fa --image_type png {params.chmfolder} {params.bins} {params.plotsfolder} {params.tetra} 95 2>> {log}
        checkm nx_plot -x fa --image_type pdf {params.bins} {params.plotsfolder} 2>> {log}
        #checkm nx_plot -x fa --image_type png {params.bins} {params.plotsfolder} 2>> {log}
        checkm len_plot -x fa --image_type pdf {params.bins} {params.plotsfolder} 2>> {log}
        #checkm len_plot -x fa --image_type png {params.bins} {params.plotsfolder} 2>> {log}
        checkm len_hist -x fa --image_type pdf {params.bins} {params.plotsfolder} 2>> {log}
        #checkm len_hist -x fa --image_type png {params.bins} {params.plotsfolder} 2>> {log}
        checkm marker_plot -x fa --image_type pdf {params.chmfolder} {params.bins} {params.plotsfolder} 2>> {log}
        #checkm marker_plot -x fa --image_type png {params.chmfolder} {params.bins} {params.plotsfolder} 2>> {log}
        checkm coverage -x fa -t {threads} {params.bins} {params.coverage} {input.bam} 2>> {log}
        checkm par_plot -x fa --image_type pdf {params.chmfolder} {params.bins} {params.plotsfolder} {params.coverage} 2>> {log}
        #checkm par_plot -x fa --image_type png {params.chmfolder} {params.bins} {params.plotsfolder} {params.coverage} 2>> {log}
        checkm profile --tab_table {params.coverage} > {params.coverage_profile} 2>> {log}       
        '''



rule eukrep:
    input:
        contigs = ASSEMBLYFOLDER + '{sample}/{sample}_spades/contigs_min5K.fasta'
    output:
        euk = EUKREPFOLDER + '{sample}/{sample}_euk_contigs_min5K.fasta',
        prok = EUKREPFOLDER + '{sample}/{sample}_prok_contigs_min5K.fasta'
    log:
        EUKREPFOLDER + '{sample}/{sample}.eukrep.log'
    benchmark:
        EUKREPFOLDER + '{sample}/{sample}.eukrep.benchmark'
    threads:
        1
    shell:
        '''
        EukRep -i {input.contigs} -o {output.euk} --prokarya {output.prok} 2> {log}
        '''

rule copyContigs:
    input:
        contigs = ASSEMBLYFOLDER + '{sample}/{sample}_spades/contigs_min5K.fasta'
    output:
        contigs = EUKREPFOLDER + '{sample}/{sample}_all_contigs_min5K.fasta'
    threads:
        1
    shell:
        '''
        cp {input.contigs} {output.contigs}
        '''

rule prodigal:
    input:
        contigs = ASSEMBLYFOLDER + '{sample}/{sample}_spades/contigs_min5K.fasta'
    output:
        marker = touch(ANNOFOLDER + '{sample}/{sample}_prodigal.done'),
        fna = ANNOFOLDER + '{sample}/{sample}.fna',
        faa = ANNOFOLDER + '{sample}/{sample}.faa',
        gff = ANNOFOLDER + '{sample}/{sample}.gff'
    log:
        ANNOFOLDER + '{sample}/{sample}_prodigal.log'
    threads:
        1
    benchmark:
        ANNOFOLDER + '{sample}/{sample}_prodigal.benchmark'
    shell:
        '''
        prodigal -a {output.faa} -d {output.fna} -f gff -o {output.gff} -c -q -m  -p meta -i {input.contigs}
        '''

rule checkm_ssu_finder:
    input:
        contigs = EUKREPFOLDER + '{sample}/{sample}_{euk_prok}_contigs_min5K.fasta',
        checkm_done = CHECKMFOLDER + '{sample}/{euk_prok}/{sample}_checkm.done'
    output:
        marker = touch(CHECKMFOLDER + '{sample}/{euk_prok}/{sample}_checkm_ssu_finder.done')
    threads:
        16
    params:
        metabatbinfolder = METABATFOLDER + '{sample}/{euk_prok}/bins/',
        ssufolder = CHECKMFOLDER + '{sample}/{euk_prok}/ssu_finder/'
    log:
        CHECKMFOLDER + '{sample}/{euk_prok}/{sample}_checkm_ssu_finder.log'
    benchmark:
        CHECKMFOLDER + '{sample}/{euk_prok}/{sample}_checkm_ssu_finder.benchmark'
    shell:
        '''
        checkm ssu_finder -x fa -t {threads} {input.contigs}  {params.metabatbinfolder} {params.ssufolder} &> {log}
        '''

