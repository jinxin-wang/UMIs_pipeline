## When you have an unmapped BAM file with RX tags, use a combination of an
## aligner and Picard’s MergeBamAlignment tool to generate a mapped BAM that
## includes all necessary metadata. 

rule merge_align:
    input:
        "bam/{sample}.unmapped.withUMI.bam"
    output:
        temp("bam/{sample}.mapped.bam")
    log:
        "logs/bam/{sample}.mapped.bam.log"
    params:
        queue = "shortq",
        index = config["REF"]["INDEX"]
    threads: 16
    resources:
        mem_mb = 51200
    conda: "UMIs"
    shell:
        "picard SamToFastq I={input} F=/dev/stdout INTERLEAVE=true "
        " | bwa mem -p -t 16 {params.index} /dev/stdin "
        " | picard MergeBamAlignment ALIGNED=/dev/stdin UNMAPPED={input} O={output} R={params.index} "
        "     SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 "
        "     ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true 2> {log} "
