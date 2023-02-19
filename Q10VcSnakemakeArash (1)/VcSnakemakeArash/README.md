# Varriant calling pipeline using Snakemake

## what is Variant calling?

Variant calling is the process of identifying the differences between a reference genome and a sample genome. The differences are called variants and can be SNPs, indels, or structural variants. Variant calling is a critical step in the analysis of next-generation sequencing data. It is the first step in the analysis of most genetic diseases and is used to identify the genetic variants that are associated with a disease. It is also used to identify the genetic variants that are associated with drug response and drug resistance.

## What is Snakemake?

Snakemake is a workflow management system for computational pipelines. It is a Python-based tool that allows for easy creation, management, and execution of complex workflows. Snakemake provides a way to define rules that specify how to produce output files from input files. These rules are then combined into a workflow graph and executed by Snakemake, which determines the order in which to execute the rules based on the dependencies between files. Snakemake is particularly useful in bioinformatics, where many tasks involve processing large amounts of data through multiple steps, each of which may depend on the output of previous steps.

## This workflow

This is a Snakemake workflow written in Python that processes NGS (next-generation sequencing) data to call variants in a set of samples.

The workflow starts by reading a tab-separated file with two columns 'sample_name' and 'fastq', which correspond to the names of the samples and the path to the fastq files, respectively. The read file is converted to a Pandas DataFrame and is set to be the index of the DataFrame with the sample_name column.

The workflow consists of seven rules:

1. all: This rule specifies the final output of the workflow, which is a list of variant call format (VCF) files for each sample.
2. bowtie_index: This rule builds the bowtie2 index for the reference genome specified in the config file 'config/config.yaml'.
3. map_reads: This rule maps the reads to the reference genome using the bowtie2 program. The output is a SAM (sequence alignment/map) file.
4. sam_to_bam: This rule converts the SAM file to a BAM (binary version of SAM) file and sorts it using the samtools program.
5. index_bam: This rule indexes the BAM file using samtools.
6. call_variants: This rule calls variants using bcftools by performing a pileup of the aligned reads and calling variants from the pileup.

All rules are executed within a conda environment specified in a separate YAML file.

## Running the workflow

1. Copy your genome fasta file into resources/genome/ directory and update the config/config.yaml genome path
2. Copy fastq read files into resources/ directory and update config/sample.tsv (tab-seperated)
3. Run the command bellow

```console

$ snakemake -p --use-conda -j2

```

4. the .vcf files can be filnd in results/01_called_variants directory

example output:


├── config

│   ├── config.yaml
│   └── samples.tsv
├── resources
│   ├── genome
│   │   ├── yeastGenome2.fa
│   │   ├── yeastGenome.fa
│   │   ├── yeastGenome.fa.1.bt2
│   │   ├── yeastGenome.fa.2.bt2
│   │   ├── yeastGenome.fa.3.bt2
│   │   ├── yeastGenome.fa.4.bt2
│   │   ├── yeastGenome.fa.fai
│   │   ├── yeastGenome.fa.rev.1.bt2
│   │   └── yeastGenome.fa.rev.2.bt2
│   └── yeastReads.fastq
├── results
│   ├── 00_mapped_reads
│   │   ├── sample_1.bam
│   │   └── sample_1.bam.bai
│   └── 01_called_variants
│       └── sample_1.vcf
└── workflow
    ├── envs
    │   ├── bowtie2.yaml
    │   └── htslib.yaml
    └── Snakefile
