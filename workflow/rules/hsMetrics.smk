## A rule to get mapping stats, with GATK CollectHsMetrics
rule collect_Hs_metrics:
    input:
        bam = "bam/{sample}.mapped.bam",
        bai = "bam/{sample}.mapped.bai",
    output:
        "HsMetrics/{sample}_HsMetrics.tsv"
    params:
        queue    = "mediumq",
        gatk     = config["gatk"]["app"],
        bait     = config["gatk"][config["samples"]][config["seq_type"]]["hsmetrics_bait"],
        interval = config["gatk"][config["samples"]][config["seq_type"]]["hsmetrics_interval"],
    log:
        "logs/mapping_QC/{sample}_HsMetrics.log"
    threads : 4
    resources:
        mem_mb = 40960
    shell:
        "{params.gatk} --java-options \"-Xmx40g  -XX:+UseParallelGC -XX:ParallelGCThreads={threads} -Djava.io.tmpdir=/mnt/beegfs/userdata/$USER/tmp \" CollectHsMetrics"
        " -TI {params.interval}"
        " -BI {params.bait}"
        " -I {input.bam}"
        " -O {output} 2> {log}"
