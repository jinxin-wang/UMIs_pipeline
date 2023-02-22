## Use fgbio’s ClipBam to eliminate overlap between reads to ensure that
## downstream processes, particularly variant calling, cannot double count evidence
## from the same template when both reads span a variant site in the same
## template. Clipping overlapping reads is only performed on “FR” read pairs and
## is implemented by clipping approximately half the overlapping bases from each
## read. By default, hard clipping is performed; soft-clipping may be substituted
## using the “--soft-clip” parameter.

rule clip_fltered_bam:
    input:
        "bam/{sample}.consensus.mapped.filtered.bam"
    output:
        "bam/{sample}.consensus.mapped.filtered.clipped.bam"
    log:
        "logs/bam/{sample}.consensus.mapped.filtered.clipped.bam.log"
    params:
        queue = "mediumq",
        index = config["REF"]["INDEX"]
    threads: 16
    resources:
        mem_mb = 51200
    conda: "UMIs"
    shell:
        "fgbio ClipBam "
        "  --input={input} "
        "  --output={output} "
        "  --ref={params.index} "
        "  --soft-clip=false "
        "  --clip-overlapping-reads=true 2> {log} "


## Variant calling can be accomplished with the variant caller of your choice. The
## following rule shows is in tumor-only mode to generate
## a VCF file, and Picard’s SortVcf to sort and index the resulting VCF.

rule produce_vcf_tumor_only_step1:
    input:
        "bam/{sample}.consensus.mapped.filtered.clipped.bam"
    output:
        temp("vcf/{sample}.tmp.vcf")
    log:
        "logs/vcf/{sample}.tmp.vcf.log"
    params:
        queue = "mediumq",
        min_af= 0.01,
        index = config["REF"]["INDEX"],
        target_region = config["REF"]["TARGET_REGION"]
    threads: 16
    resources:
        mem_mb = 51200
    conda: "UMIs"
    shell:
        "vardict-java "
        "  -G {params.index} "
        "  -N {sample} "
        "  -f {params.min_af} "
        "  -b {input} "
        "  -z -c 1 -s 2 -E 3 -g 4 -th 4 "
        "  {params.target_region} "
        "  | teststrandbias.R "
        "  | var2vcf_valid.pl -N {sample} -E -f {params.min_af} "
        "  awk '{if ($1 ~ /^#/) print; else if ($4 != $5) print}' " 
        "  > {output} 2> {log}"

rule produce_vcf_tumor_only_step2:
    input:
        "vcf/{sample}.tmp.vcf"
    output:
        "vcf/{sample}.vcf"
    log:
        "logs/bam/{sample}.vcf.log"
    params:
        queue = "mediumq",
        min_af= 0.01,
        dict = config["REF"]["DICT"],
    threads: 16
    resources:
        mem_mb = 51200
    conda: "UMIs"
    shell:
        " picard SortVcf "
        "  I={input} "
        "  O={output} "
        "  SD={params.dict} 2> {log} "
