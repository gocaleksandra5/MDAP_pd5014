configfile: "config.yaml"

SAMPLE = config["sample"]
R1 = config["reads"]["R1"]
R2 = config["reads"]["R2"]
REF = config["reference"]
ADAPTER = config["adapter"]
OUT = config["output"]
THREADS = config["threads"]


rule all:
    input:
        f"{OUT}/fastqc_raw/{SAMPLE}_1_fastqc.html",
        f"{OUT}/fastqc_raw/{SAMPLE}_2_fastqc.html",

        f"{OUT}/trimmed/{SAMPLE}_1_paired.fastq.gz",
        f"{OUT}/trimmed/{SAMPLE}_2_paired.fastq.gz",

        f"{OUT}/fastqc_trimmed/{SAMPLE}_1_paired_fastqc.html",
        f"{OUT}/fastqc_trimmed/{SAMPLE}_2_paired_fastqc.html",

        f"{OUT}/bam/{SAMPLE}_sorted.bam",
        f"{OUT}/bam/{SAMPLE}_sorted.bam.bai",

        f"{OUT}/coverage/{SAMPLE}_coverage.txt",

        f"{OUT}/variants/{SAMPLE}.vcf",

        f"{OUT}/multiqc/multiqc_report.html"


rule fastqc_raw:
    input:
        R1,
        R2
    output:
        html1 = f"{OUT}/fastqc_raw/{SAMPLE}_1_fastqc.html",
        zip1 = f"{OUT}/fastqc_raw/{SAMPLE}_1_fastqc.zip",
        html2 = f"{OUT}/fastqc_raw/{SAMPLE}_2_fastqc.html",
        zip2 = f"{OUT}/fastqc_raw/{SAMPLE}_2_fastqc.zip"
    threads:
        THREADS
    shell:
        """
        mkdir -p {OUT}/fastqc_raw

        fastqc \
            -t {threads} \
            -o {OUT}/fastqc_raw \
            {input}
        """


rule trim:
    input:
        R1 = R1,
        R2 = R2,
        adapter = ADAPTER
    output:
        R1_paired = f"{OUT}/trimmed/{SAMPLE}_1_paired.fastq.gz",
        R1_unpaired = f"{OUT}/trimmed/{SAMPLE}_1_unpaired.fastq.gz",
        R2_paired = f"{OUT}/trimmed/{SAMPLE}_2_paired.fastq.gz",
        R2_unpaired = f"{OUT}/trimmed/{SAMPLE}_2_unpaired.fastq.gz"
    threads:
        THREADS
    params:
        leading = config["trimmomatic"]["leading"],
        trailing = config["trimmomatic"]["trailing"],
        slidingwindow = config["trimmomatic"]["slidingwindow"],
        minlen = config["trimmomatic"]["minlen"]
    shell:
        """
        mkdir -p {OUT}/trimmed

        trimmomatic PE \
            -threads {threads} \
            {input.R1} \
            {input.R2} \
            {output.R1_paired} \
            {output.R1_unpaired} \
            {output.R2_paired} \
            {output.R2_unpaired} \
            ILLUMINACLIP:{input.adapter}:2:30:10 \
            LEADING:{params.leading} \
            TRAILING:{params.trailing} \
            SLIDINGWINDOW:{params.slidingwindow} \
            MINLEN:{params.minlen}
        """


rule fastqc_trimmed:
    input:
        R1 = f"{OUT}/trimmed/{SAMPLE}_1_paired.fastq.gz",
        R2 = f"{OUT}/trimmed/{SAMPLE}_2_paired.fastq.gz"
    output:
        html1 = f"{OUT}/fastqc_trimmed/{SAMPLE}_1_paired_fastqc.html",
        zip1 = f"{OUT}/fastqc_trimmed/{SAMPLE}_1_paired_fastqc.zip",
        html2 = f"{OUT}/fastqc_trimmed/{SAMPLE}_2_paired_fastqc.html",
        zip2 = f"{OUT}/fastqc_trimmed/{SAMPLE}_2_paired_fastqc.zip"
    threads:
        THREADS
    shell:
        """
        mkdir -p {OUT}/fastqc_trimmed

        fastqc \
            -t {threads} \
            -o {OUT}/fastqc_trimmed \
            {input}
        """


rule bwa_index:
    input:
        REF
    output:
        amb = REF + ".amb",
        ann = REF + ".ann",
        bwt = REF + ".bwt",
        pac = REF + ".pac",
        sa = REF + ".sa"
    shell:
        """
        bwa index {input}
        """

rule faidx:
    input:
        REF
    output:
        REF + ".fai"
    shell:
        """
        samtools faidx {input}
	"""

rule mapping:
    input:
        ref = REF,
        index = rules.bwa_index.output,
        R1 = f"{OUT}/trimmed/{SAMPLE}_1_paired.fastq.gz",
        R2 = f"{OUT}/trimmed/{SAMPLE}_2_paired.fastq.gz"
    output:
        temp(f"{OUT}/bam/{SAMPLE}.bam")
    threads:
        THREADS
    shell:
        """
        mkdir -p {OUT}/bam

        bwa mem \
            -t {threads} \
            {input.ref} \
            {input.R1} \
            {input.R2} \
        | samtools view \
            -b \
            -o {output}
        """


rule sort_bam:
    input:
        f"{OUT}/bam/{SAMPLE}.bam"
    output:
        f"{OUT}/bam/{SAMPLE}_sorted.bam"
    threads:
        THREADS
    shell:
        """
        samtools sort \
            -@ {threads} \
            -o {output} \
            {input}
        """


rule index_bam:
    input:
        f"{OUT}/bam/{SAMPLE}_sorted.bam"
    output:
        f"{OUT}/bam/{SAMPLE}_sorted.bam.bai"
    shell:
        """
        samtools index {input}
        """


rule coverage:
    input:
        bam = f"{OUT}/bam/{SAMPLE}_sorted.bam",
        bai = f"{OUT}/bam/{SAMPLE}_sorted.bam.bai"
    output:
        f"{OUT}/coverage/{SAMPLE}_coverage.txt"
    shell:
        """
        mkdir -p {OUT}/coverage

        samtools depth \
            -a \
            {input.bam} \
            > {output}
        """


rule variants:
    input:
        ref = REF,
        bam = f"{OUT}/bam/{SAMPLE}_sorted.bam",
        bai = f"{OUT}/bam/{SAMPLE}_sorted.bam.bai"
    output:
        f"{OUT}/variants/{SAMPLE}.vcf"
    shell:
        """
        mkdir -p {OUT}/variants

        bcftools mpileup \
            -Ou \
            -f {input.ref} \
            {input.bam} \
        | bcftools call \
            -mv \
            -Ov \
            -o {output}
        """


rule multiqc:
    input:
        rules.fastqc_raw.output,
        rules.fastqc_trimmed.output
    output:
        f"{OUT}/multiqc/multiqc_report.html"
    shell:
        """
        mkdir -p {OUT}/multiqc

        multiqc \
            {OUT}/fastqc_raw \
            {OUT}/fastqc_trimmed \
            -o {OUT}/multiqc \
            -n multiqc_report.html \
            -f
        """
