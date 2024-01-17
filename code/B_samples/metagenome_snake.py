# To run:
# snakemake -k -s code/metagenome_snake.py --config norm=False --cluster "DIR=\$(dirname {params.errFile}); mkdir -p \"\${{DIR}}\"; qsub -cwd -S /bin/bash -N snakejob -V -pe {params.pe} -l h_vmem={params.mem} -e {params.errFile} -o {params.outFile}" --jobs 8

import os

samples = ["20200831.B-CP_Fresh_Pro", "20200831.B-MkB_Fresh_Pro", "20200831.B-MK_Sed", "20200831.B-OP_Sed", "20200831.B-CP_Sed", "20200831.B-MkP_Fresh_Pro", "20200831.B-OP_Fresh_Pro"]
gtdbtks = ["scratch/{}/annotation/{}.bac120.classify.tree".format(sample,sample) for sample in samples]
checkms = ["scratch/{}/checkm/{}.checkm.done".format(sample,sample) for sample in samples]
eukreps = ["scratch/{}/eukrep/pro_contigs.fasta".format(sample) for sample in samples]
smashes = ["scratch/{}/antismash/results/index.html".format(sample) for sample in samples]
ftables = ["scratch/{}/antismash/bin_features_table.txt".format(sample) for sample in samples]

rule all:
    input:
        gtdbtks,
        checkms,
        eukreps,
        smashes,
        ftables

rule qc:
    input:
        r1 = "data/{sample}.1.fq.gz",
        r2 = "data/{sample}.2.fq.gz"
    output:
        qc1 = "scratch/{sample}/qc/{sample}.1.fq.gz",
        qc2 = "scratch/{sample}/qc/{sample}.2.fq.gz",
        fail = "scratch/{sample}/qc/{sample}_fail.fq.gz",
        single = "scratch/{sample}/qc/{sample}_single.fq.gz",
        stats = "scratch/{sample}/qc/{sample}_qc.stats",
        adaptMatch = "scratch/{sample}/qc/{sample}_adaptMatch.fq.gz",
        adaptSingle = "scratch/{sample}/qc/{sample}_adaptSingle.fq.gz",
        adaptStats = "scratch/{sample}/qc/{sample}_adapt.stats",
        phixMatch = "scratch/{sample}/qc/{sample}_phixMatch.fq.gz",
        phixSingle = "scratch/{sample}/qc/{sample}_phixSingle.fq.gz",
        phixStats = "scratch/{sample}/qc/{sample}_phix.stats"
    params:
        errFile = "scratch/{sample}/qc/qc.err.log",
        outFile = "scratch/{sample}/qc/qc.out.log",
        pe = "smp 32",
        mem = "2G"
    log:
        log = "scratch/{sample}/qc/{sample}.qc.log",
        adaptLog = "scratch/{sample}/qc/{sample}.adapt.log",
        phixLog = "scratch/{sample}/qc/{sample}.phix.log"
    shell:
        """
        ml BBMap
        bbduk.sh -Xmx2G usejni=t threads=14 overwrite=t qin=33 in1={input.r1} in2={input.r2} ref=/nfs/modules/modules/software/BBMap/38.26-foss-2018b/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 out=stdout.fq outm={output.adaptMatch} outs={output.adaptSingle} refstats={output.adaptStats} statscolumns=5 2> {log.adaptLog} | bbduk.sh -Xmx2G usejni=t threads=4 interleaved=true overwrite=t qin=33 in=stdin.fq out=stdout.fq outm={output.phixMatch} outs={output.phixSingle} ref=/nfs/modules/modules/software/BBMap/38.26-foss-2018b/resources/phix174_ill.ref.fa.gz k=31 hdist=1 refstats={output.phixStats} statscolumns=5 2> {log.phixLog} | bbduk.sh -Xmx2G usejni=t threads=14 overwrite=t interleaved=true qin=33 in=stdin.fq fastawrap=10000 out1={output.qc1} out2={output.qc2} interleaved=true outm={output.fail} outs={output.single} minlength=45 maq=20 maxns=0 overwrite=t stats={output.stats} statscolumns=5 trimq=14 qtrim=r 2>{log.log}
        ml purge
        """

if config["norm"]:
    rule norm:
        input:
            qc1 = "scratch/{sample}/qc/{sample}.1.fq.gz",
            qc2 = "scratch/{sample}/qc/{sample}.2.fq.gz",
            single = "scratch/{sample}/qc/{sample}_single.fq.gz"
        output:
            n1 = "scratch/{sample}/norm/{sample}.1.fq.gz",
            n2 = "scratch/{sample}/norm/{sample}.2.fq.gz",
            ns = "scratch/{sample}/norm/{sample}_single.fq.gz",
            histPEin = "scratch/{sample}/norm/{sample}_pe.hist",
            histSin = "scratch/{sample}/norm/{sample}_s.hist",
            histPEout = "scratch/{sample}/norm/{sample}_pe.norm.hist",
            histSout = "scratch/{sample}/norm/{sample}_s.norm.hist"
        params:
            errFile = "scratch/{sample}/norm/norm.err.log",
            outFile = "scratch/{sample}/norm/norm.out.log",
            pe = "smp 32",
            mem = "8G"
        log:
            log = "scratch/{sample}/norm/{sample}.norm.log"
        shell:
            """
            ml BBMap
            bbnorm.sh -Xmx200G in={input.qc1} in2={input.qc2} extra={input.single} out={output.n1} out2={output.n2} target=80 mindepth=1 hist={output.histPEin} histout={output.histPEout}
            bbnorm.sh -Xmx200G extra={input.qc1} extra={input.qc2} in={input.single} out={output.ns} target=80 mindepth=1 hist={output.histSin} histout={output.histSout}
            ml purge
            """

    rule assembly:
        input:
            n1 = "scratch/{sample}/norm/{sample}.1.fq.gz",
            n2 = "scratch/{sample}/norm/{sample}.2.fq.gz",
            ns = "scratch/{sample}/norm/{sample}_single.fq.gz"
        output:
            done = touch("scratch/{sample}/assembly/{sample}.ass.done")
        params:
            outdir = "scratch/{sample}/assembly",
            errFile = "scratch/{sample}/assembly/spades.err.log",
            outFile = "scratch/{sample}/assembly/spades.out.log",
            pe = "smp 32",
            mem = "16G"
        log:
            log = "scratch/{sample}/assembly/{sample}.ass.log"
        shell:
            """
            ml SPAdes
            spades.py --meta --only-assembler --pe1-1 {input.n1} --pe1-2 {input.n2} --pe1-s {input.ns} -o {params.outdir} -t 32 -m 512
            ml purge
            """
else:
    rule assembly:
        input:
            qc1 = "scratch/{sample}/qc/{sample}.1.fq.gz",
            qc2 = "scratch/{sample}/qc/{sample}.2.fq.gz",
            single = "scratch/{sample}/qc/{sample}_single.fq.gz"
        output:
            done = touch("scratch/{sample}/assembly/{sample}.ass.done"),
            scaffolds = "scratch/{sample}/assembly/contigs.fasta"
        params:
            outdir = "scratch/{sample}/assembly",
            errFile = "scratch/{sample}/assembly/spades.err.log",
            outFile = "scratch/{sample}/assembly/spades.out.log",
            pe = "smp 64",
            mem = "48G"
        log:
            log = "scratch/{sample}/assembly/{sample}.ass.log"
        shell:
            """
            ml SPAdes
            spades.py --meta --only-assembler --pe1-1 {input.qc1} --pe1-2 {input.qc2} --pe1-s {input.single} -o {params.outdir} -t 64 -m 3072
            ml purge
            """

rule filter:
    input:
        done = "scratch/{sample}/assembly/{sample}.ass.done",
        contigs = "scratch/{sample}/assembly/contigs.fasta"
    output:
        contigs = "scratch/{sample}/assembly/{sample}_contigs.1k.fasta"
    params:
        errFile = "scratch/{sample}/assembly/filter.err.log",
        outFile = "scratch/{sample}/assembly/filter.out.log",
        pe = "smp 1",
        mem = "8G"
    shell:
        """
        ml Perl
        perl /nfs/home/fieldc/scripts/contig-filter.pl -f {input.contigs} -c 1000 > {output.contigs}
        ml purge
        """

rule align:
    input:
        qc1 = "scratch/{sample}/qc/{sample}.1.fq.gz",
        qc2 = "scratch/{sample}/qc/{sample}.2.fq.gz",
        single = "scratch/{sample}/qc/{sample}_single.fq.gz",
        contigs = "scratch/{sample}/assembly/{sample}_contigs.1k.fasta"
    output:
        sampe = temp("scratch/{sample}/metabat/{sample}_contigs_pe.1k.sam"),
        bampe = temp("scratch/{sample}/metabat/{sample}_contigs_pe.1k.bam"),
        sams = temp("scratch/{sample}/metabat/{sample}_contigs_s.1k.sam"),
        bams = temp("scratch/{sample}/metabat/{sample}_contigs_s.1k.bam"),
        bam = "scratch/{sample}/metabat/{sample}_contigs.1k.bam",
        done = touch("scratch/{sample}/metabat/{sample}.align.done")
    params:
        errFile = "scratch/{sample}/metabat/align.err.log",
        outFile = "scratch/{sample}/metabat/align.out.log",
        pe = "smp 32",
        mem = "8G"
    shell:
        """
        ml BWA
        ml SAMtools
        bwa index {input.contigs}
        bwa mem -t 32 {input.contigs} {input.qc1} {input.qc2} | samtools view -h -F 4 - | python /nfs/cds/mOTUv2/mOTUs_v2/bin/msamtools_python.py 0.97 200 90 None > {output.sampe}
        bwa mem -t 32 {input.contigs} {input.single} | samtools view -h -F 4 - | python /nfs/cds/mOTUv2/mOTUs_v2/bin/msamtools_python.py 0.97 200 90 None > {output.sams}
        samtools view -u {output.sampe} | samtools sort -@ 32 - > {output.bampe}
        samtools view -u {output.sams} | samtools sort -@ 32 - > {output.bams}
        samtools merge -u -@ 32 - {output.bampe} {output.bams} | samtools sort -@ 32 - > {output.bam}
        samtools index {output.bam}
        ml purge
        """

rule metabat:
    input:
        done = "scratch/{sample}/metabat/{sample}.align.done",
        contigs = "scratch/{sample}/assembly/{sample}_contigs.1k.fasta",
        bam = "scratch/{sample}/metabat/{sample}_contigs.1k.bam"
    output:
        depth = "scratch/{sample}/metabat/{sample}.depth",
        done = touch("scratch/{sample}/metabat/{sample}.metabat.done")
    params:
        outdir = "scratch/{sample}/metabat/{sample}",
        errFile = "scratch/{sample}/metabat/metabat.err.log",
        outFile = "scratch/{sample}/metabat/metabat.out.log",
        pe = "smp 32",
        mem = "8G"
    shell:
        """
        ml MetaBAT
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}
        metabat -t 32 -i {input.contigs} -o {params.outdir} --unbinned -v -a {output.depth}
        ml purge
        """

rule checkm:
    input:
        done = "scratch/{sample}/metabat/{sample}.metabat.done",
        contigs = "scratch/{sample}/assembly/{sample}_contigs.1k.fasta",
        bam = "scratch/{sample}/metabat/{sample}_contigs.1k.bam"
    output:
        tetra = "scratch/{sample}/checkm/{sample}.tetra.tsv",
        cov = "scratch/{sample}/checkm/{sample}.coverage.tsv",
        profile = "scratch/{sample}/checkm/{sample}.profile.txt",
        table = "scratch/{sample}/checkm/{sample}.results.txt",
        done = touch("scratch/{sample}/checkm/{sample}.checkm.done")
    params:
        metabatdir = "scratch/{sample}/metabat/",
        outputdir = "scratch/{sample}/checkm/output",
        plotsdir = "scratch/{sample}/checkm/plots",
        errFile = "scratch/{sample}/checkm/checkm.err.log",
        outFile = "scratch/{sample}/checkm/checkm.out.log",
        pe = "smp 1",
        mem = "128G"
    shell:
        """
        ml CheckM
        checkm lineage_wf -t 1 -x fa -f {output.table} --tab_table {params.metabatdir} {params.outputdir}
        checkm tetra {input.contigs} {output.tetra}
        checkm bin_qa_plot -x fa {params.outputdir} {params.metabatdir} {params.plotsdir}
        checkm dist_plot -x fa {params.outputdir} {params.metabatdir} {params.plotsdir} {output.tetra} 95
        checkm nx_plot -x fa {params.metabatdir} {params.plotsdir}
        checkm len_plot -x fa {params.metabatdir} {params.plotsdir}
        checkm len_hist -x fa {params.metabatdir} {params.plotsdir}
        #checkm marker_plot -x fa {params.outputdir} {params.metabatdir} {params.plotsdir} --dpi 300
        checkm coverage -x fa -t 1 {params.metabatdir} {output.cov} {input.bam}
        checkm par_plot -x fa {params.outputdir} {params.metabatdir} {params.plotsdir} {output.cov}
        checkm profile --tab_table {output.cov} > {output.profile}
        ml purge
        """

rule prodigal:
    input:
        contigs = "scratch/{sample}/assembly/{sample}_contigs.1k.fasta"
    output:
        fna = "scratch/{sample}/annotation/{sample}.fna",
        faa = "scratch/{sample}/annotation/{sample}.faa",
        gff = "scratch/{sample}/annotation/{sample}.gff"
    params:
        errFile = "scratch/{sample}/annotation/prodigal.err.log",
        outFile = "scratch/{sample}/annotation/prodigal.out.log",
        pe = "smp 16",
        mem = "4G"
    shell:
        """
        ml prodigal
        prodigal -a {output.faa} -d {output.fna} -f gff -o {output.gff} -c -q -m -p meta -i {input.contigs}
        ml purge
        """

rule gtdbtk:
    input:
        done = "scratch/{sample}/metabat/{sample}.metabat.done"
    output:
        tree = "scratch/{sample}/annotation/{sample}.bac120.classify.tree"
    params:
        metabatdir = "scratch/{sample}/metabat/",
        outputdir = "scratch/{sample}/annotation/",
        errFile = "scratch/{sample}/annotation/gtdbtk.err.log",
        outFile = "scratch/{sample}/annotation/gtdbtk.out.log",
        prefix = "{sample}",
        pe = "smp 1",
        mem = "512G"
    shell:
        """
        ml GTDBTk
        gtdbtk classify_wf --genome_dir {params.metabatdir} -x fa --out_dir {params.outputdir} --cpus 1 --prefix {params.prefix}
        ml purge
        """

rule eukrep:
    input:
        contigs = "scratch/{sample}/assembly/{sample}_contigs.1k.fasta"
    output:
        euk_contigs = "scratch/{sample}/eukrep/euk_contigs.fasta",
        pro_contigs = "scratch/{sample}/eukrep/pro_contigs.fasta"
    params:
        pe = "smp 8",
        mem = "8G",
        errFile = "scratch/{sample}/eukrep/eukrep.err.log",
        outFile = "scratch/{sample}/eukrep/eukrep.out.log"
    shell:
        """
        ml EukRep
        EukRep -i {input.contigs} -o {output.euk_contigs} --prokarya {output.pro_contigs} --seq_names
        ml purge
        """

rule antismash:
    input:
        contigs = "scratch/{sample}/assembly/{sample}_contigs.1k.fasta",
        gff = "scratch/{sample}/annotation/{sample}.gff"
    output:
        webpage = "scratch/{sample}/antismash/results/index.html"
    params:
        outputdir = "scratch/{sample}/antismash/results",
        pe = "smp 32",
        mem = "8G",
        errFile = "scratch/{sample}/antismash/antismash.err.log",
        outFile = "scratch/{sample}/antismash/antismash.out.log"
    shell:
        """
        ml antiSMASH
        antismash -c 32 --output-dir {params.outputdir} --genefinding-gff3 {input.gff} {input.contigs} --minlength 5000
        ml purge
        """

rule postprocessing:
    input:
        webpage = "scratch/{sample}/antismash/results/index.html",
    output:
        feature_table = "scratch/{sample}/antismash/bin_features_table.txt",
    params:
        dir = "{sample}",
        pe = "smp 1",
        mem = "8G",
        errFile = "scratch/{sample}/antismash/postprocessing.err.log",
        outFile = "scratch/{sample}/antismash/postprocessing.out.log",
    shell:
        """
        ./code/metagenome_postprocessing.sh {params.dir}
        """

