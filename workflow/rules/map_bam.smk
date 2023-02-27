## 1. Identify which reads come from the same source molecule by using fgbio’s
##    GroupReadsByUmi tool, which assigns a unique source molecule ID to each
##    applicable read, stores the ID in the MI tag, and outputs a BAM file that is sorted
##    by the MI tag and ready for consensus calling. The source molecule is identified
##    using a combination of UMI sequence and mapping positions from reads 1
##    and 2.

##    GroupReadsByUmi implements several strategies for matching UMIs to account
##    for sequencing error. The paired method implements the directed adjacency
##    graph method introduced by UMI-tools [1]. Parameters are available to control
##    how many errors are allowed when matching UMIs at the same position and for
##    filtering (i.e., ignoring) reads with low mapping quality. Reads with low mapping
##    quality should be ignored to prevent multiple consensus reads from being
##    generated from multiple mismapped copies of the same source molecule.

rule group_reads:
    input:
        "bam/{sample}.mapped.bam"
    output:
        temp("bam/{sample}.grouped.bam")
    log:
        "logs/bam/{sample}.grouped.bam.log"
    params:
        queue = "shortq",
    threads: 16
    resources:
        mem_mb = 51200
    conda: "UMIs"
    shell:
        "fgbio GroupReadsByUmi "
        "    --input={input} "
        "    --output={output} "
        "    --strategy=paired "
        "    --edits=1 "
        "    --min-map-q=20 2> {log} "


## 2. Combine each set of reads to generate consensus reads using fgbio’s
##    CallDuplexConsensusReads. This step generates unmapped consensus reads
##    from the output of GroupReadsByUmi. There are many parameters that affect
##    the consensus calling; for an up-to-date listing and supporting documentation,
##    run CallDuplexConsensusReads with the -h option.

##    This script produces consensus reads for all molecules that have at least one
##    observation.

##    Note: CallDuplexConsensusReads requires the input BAM file to be in
##    "TemplateCoordinate" sort order. BAM file outputs by GroupReadsByUmi will
##    have the correct sort order. If needed, use the SortSam function in fgbio with
##    the argument "--sort-order=TemplateCoordinate" to format your BAM file.

rule combine_reads:
    input:
        "bam/{sample}.grouped.bam"
    output:
        temp("bam/{sample}.consensus.unmapped.bam")
    log:
        "logs/bam/{sample}.consensus.unmapped.bam.log"
    params:
        queue = "shortq",
    threads: 16
    resources:
        mem_mb = 51200
    conda: "UMIs"
    shell:
        "fgbio CallDuplexConsensusReads "
        "    --input={input} "
        "    --output={output} "
        "    --error-rate-pre-umi=45 "
        "    --error-rate-post-umi=30 "
        "    --min-input-base-quality=30 2> {log} "


## 3. Remap consensus reads. After you have generated consensus reads, you must
##    remap the reads. The mapping procedure is the same as for raw reads described
##    in section B:

rule remap_consensus_reads:
    input:
        "bam/{sample}.consensus.unmapped.bam"
    output:
        temp("bam/{sample}.consensus.mapped.bam")
    log:
        "logs/bam/{sample}.consensus.mapped.bam.log"
    params:
        queue = "shortq",
        index = config["REF"]["INDEX"],
        picard= config["APP"]["PICARD"],
    threads: 16
    resources:
        mem_mb = 51200
    conda: "UMIs"
    shell:
        "java -Xmx8g -jar {params.picard} SamToFastq "
        " I={input} F=/dev/stdout INTERLEAVE=true "
        " | bwa mem -p -t 16 {params.index} /dev/stdin "
        " | java -Xmx32g -jar {params.picard} MergeBamAlignment ALIGNED=/dev/stdin "
        "     UNMAPPED={input} O={output} R={params.index} "
        "     SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 "
        "     ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true 2> {log} "

## 4. Filter consensus reads using fgbio’s FilterConsensusReads. There are
##    two kinds of filtering: 1) masking or filtering individual bases in reads, and
##    2) filtering reads (i.e., not writing them to the output file). Base-level masking/
##    filtering is only applied if per-base tags are present (see the documentation for
##    CallDuplexConsensusReads for tag descriptions). Read-level filtering is always
##    applied.

##    When filtering reads, secondary alignments and supplementary records may
##    be removed independently if they fail one or more filters. If either R1 or R2
##    primary alignments fail a filter, all records for the template will be filtered out.
##    There are many parameters that affect the filtering of the consensus reads. For
##    an up-to-date listing and supporting documentation, run FilterConsensusRead
##    with the -h option. FilterConsensusRead can be applied to either mapped or
##    unmapped BAM files.

##    This script produces a filtered, consensus BAM file containing sequences
##    from molecules that have at least 3 reads for constructing the single-strand
##    consensus and have at least 6 reads for constructing the duplex consensus.

rule filter_consensus_reads:
    input:
        "bam/{sample}.consensus.mapped.bam"
    output:
        temp("bam/{sample}.consensus.mapped.filtered.bam")
    log:
        "logs/bam/{sample}.consensus.mapped.filtered.bam.log"
    params:
        queue = "shortq",
        index = config["REF"]["INDEX"]
    threads: 16
    resources:
        mem_mb = 51200
    conda: "UMIs"
    shell:
        "fgbio FilterConsensusReads "
        "  --input={input} "
        "  --output={output} "
        "  --ref={params.index} "
        "  --min-reads=6 3 3 "
        "  --max-read-error-rate=0.05 "
        "  --max-base-error-rate=0.1 "
        "  --min-base-quality=50 "
        "  --max-no-call-fraction=0.05 2> {log} "
