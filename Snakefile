#!/usr/bin/env python

configfile: "config.yaml"

SAMPLES = config["samples"]

localrules: make_barcode_file

rule all:
    input:
        "config.yaml",
        expand("fastq/{sample}.fastq.gz", sample=SAMPLES),
        expand("fastq/nontrimmed/{sample}_nontrimmed.fastq.gz", sample=SAMPLES)

# barcodes include the 'A' tail
rule make_barcode_file:
    output:
        "fastq/barcodes.tsv"
    run:
        with open(output[0], "w") as out:
            for sample, barcode in SAMPLES.items():
                out.write(f'{sample}\t{barcode}T\n')

# demultiplex with fastq-multx (10x faster than cutadapt demultiplexing)
# allow one mismatch to barcode
# remove barcodes, including 'A' tail
rule demultiplex:
    input:
        fastq = config["fastq"],
        barcodes = "fastq/barcodes.tsv"
    output:
        expand("fastq/{sample}.fastq.gz", sample=["unmatched"] + list(SAMPLES.keys())),
    log:
        "logs/demultiplex.log"
    shell: """
       (fastq-multx -B {input.barcodes} -b -m 1 {input.fastq} -o fastq/%.fastq.gz) &> {log}
        """

# check for barcode, as above, but do not trim
# these files are useful for GEO submission, which requires demultiplexed but unmodified fastq files
rule check_barcodes_only:
    input:
        fastq = config["fastq"],
        barcodes = "fastq/barcodes.tsv"
    output:
        expand("fastq/nontrimmed/{sample}_nontrimmed.fastq.gz", sample=["unmatched"] + list(SAMPLES.keys()))
    threads: config["threads"]
    log:
        "logs/check_barcodes_only.log"
    shell: """
       (fastq-multx -B {input.barcodes} -b -x -m 1 {input.fastq} -o fastq/nontrimmed/%_nontrimmed.fastq.gz) &> {log}
        """

