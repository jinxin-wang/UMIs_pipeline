## fgbio's ExtractUmisFromBam processes unmapped BAM files, where the UMIs are
## contained within read 1 or read 2 sequences, and extracts the UMI sequences into
## the RX tag.

## the first "3M2S146T" represents the the structure of one strand,
## the second "3M2S146T" represents the structure of the other strand
## "3M" represents 3 UMI bases
## "2S" represents 2 skipped bases
## "146T" represents 146 bases in the read

rule tag_with_umis:
    input:
        "bam/{sample}.unmapped.bam"
    output:
        temp("bam/{sample}.unmapped.withUMI.bam")
    log:
        "logs/bam/{sample}.unmapped.withUMI.bam.log"
    params:
        queue = "shortq",
    threads: 16
    resources:
        mem_mb = 51200
    conda: "UMIs"
    shell:
        "fgbio ExtractUmisFromBam "
        " --input={input} "
        " --output={output} "
        " --read-structure=3M2S146T 3M2S146T "
        " --molecular-index-tags=ZA ZB "
        " --single-tag=RX 2> {log} "
