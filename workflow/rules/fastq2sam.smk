## The raw data format (e.g., FASTQ or BAM) and existing processing pipelines
## determine the optimal method to generate a mapped BAM file that contains UMI
## information. Ultimately, you should extract the in-line UMI sequence from each
## read and store it in the RX tag in an unmapped BAM file. Read 1 and read 2 of the
## same pair should have the same value stored in the RX tag, regardless of which
## read contains the UMI sequence.

## If raw data are in FASTQ files instead of unmapped BAM files, use Picard’s
## FastqToSam tool to convert the FASTQ files into unmapped BAM files.

rule fastq_to_sam:
    input:
        fastq_1="DNA_samples/{sample}_R1.fastq.gz",
        fastq_2="DNA_samples/{sample}_R2.fastq.gz",
    output:
        temp("bam/{sample}.unmapped.bam")
    log:
        "logs/sam/{sample}.unmapped.bam.log"
    params:
        queue = "mediumq",
        sample= "{sample}",
    threads: 16
    resources:
        mem_mb = 51200
    conda: "UMIs"
    shell:
        "picard FastqToSam "
        " F1={input.fastq_1} "
        " F2={input.fastq_2} "
        " OUTPUT={output} "
        " SAMPLE_NAME={params.sample} > {log}"
