##### Orginal snakemake

```console

# snake make file zero
rule all:
	input:
		's1.2'

rule first:
	input:
		's1.0'
	output:
		's1.1'
	shell:
		'touch {output} '

rule second:
	input:
		's1.1'
	output:
		's1.2'
	shell:
		'touch {output}' 
```

### Q1

* Try running the above code. Do so by running snakemake in the
  example_workflow directory. What error do you get? Once you fix the
  first error, what does the second one mean?

```console
$ snakemake 
```

```
`Error: you need to specify the maximum number of CPU cores to be used at the same time. If you want to use N cores, say --cores N or -cN. For all cores on your system (be sure that this is appropriate) use --cores all. For no parallelization use --cores 1 or -c1. <_io.TextIOWrapper name='<stderr>' mode='w' encoding='utf-8'>`


```

```console
$ snakemake --cores 1
```

```
`Building DAG of jobs... MissingInputException in rule first in file /home/inf-48-2022/example_workflow/Snakefile, line 5: Missing input files for rule first:     output: s1.1     affected files:         s1.0`
```

-- the second error means that the input file is missing for the first rule and we should create it

### Q2

* Add the missing file by touching it. What happens if you run snakemake
  -n? Try re-running the workflow. What happens?

```console
$ touch s1.0
```

```console
$ snakemake -n
```

snakemake -n shows the commands that will be executed and the order in which they will be executed
so the first rule will be executed first and the second rule will be executed second.
it doesn't execute the commands, it just shows them so no files are created

```console
$ snakemake --cores 1
```

result:
it creates the file s1.1 and the file s1.2

### Q3

* Try running the workflow again. What happens? What happens if you
  remove or touch s1.1, and try re-running it?

```console
$ snakemake --cores 1
```

output:

```
`Building DAG of jobs... Nothing to be done (all requested files are present and up to date). Complete log: .snakemake/log/2023-01-29T120339.172082.snakemake.log`
```

```console
$ rm s1.1
$ snakemake --cores 1
```

the output is the same as before, nothing to be done, becuase s1.2 file is already created so snakemake dont rerun the workflow

### Q4

* Using what you’ve learned above regarding wildcards, update this workflow to be able to handle generating any file ending in .2, including the
  intermediate step as in the workflow above. Use this updated workflow
  to generate my_first_workflow.2 by running snakemake. Also use it
  to generate manual_target.2 by running snakemake manual_target.2.
  Remember that you must first generate the necessary input files manually

```console
rule all:
	input:
		'{sample}.2'

rule first:
	input:
		'{sample}.0'
	output:
		'{sample}.1'
	shell:
		'touch {output} '

rule second:
	input:
		'{sample}.1'
	output:
		'{sample}.2'
	shell:
		'touch {output}' 
```

```console
$ touch my_first_workflow.0
$ snakemake --cores 1 my_first_workflow.2
```

```console
$ touch manual_target.0
$ snakemake --cores 1 manual_target.2
```

### Q5

* Update the above workflow to use glob_wildcards or a tsv file to identify which files should be generated. Explain what effect calling expand in other rules than all has on this workflow. (Comparing dry-run outputs for different workflows with multiple targets may be helpful here.) The code solution to this questions should be handed in.

expand is used to expand the wildcards in the input and output files of the rules
so if we use expand in the all rule, it will expand the wildcards in the input files of the other rules

```console
SAMPLE, = glob_wildcards('{sample}.0')

rule all:
	input:
		expand('{sample}.2',sample=SAMPLE)

rule first:
	input:
		'{sample}.0'
	output:
		'{sample}.1'
	shell:
		'touch {output} '

rule second:
	input:
		'{sample}.1'
	output:
		'{sample}.2'
	shell:
		'touch {output}' 
```

```console
$ snakemake --cores 1
```

output:

The glob_wildcards function is used to get the list of files that match the pattern {sample}.0
and the expand function is used to expand the wildcards in the input and output files of the rules

we will have the same number of *.2 files as the number of *.0 files present in the directory

### Q6

The *-j* option is for specifying the number of cores to use for parallel processing, here 2.

The '--prioritize' or '-P' option is used in Snakemake to assign high priority to the specified targets and their dependencies in the workflow. This means that these targets and their dependencies will be executed before any other targets in the workflow.

# q7

*> multiext*:

generates multiple outputs with different file extensions. It takes two arguments: the first is the base file name (in this case, `config['genome']`), and the second is a list of extensions to add to the base file name to generate multiple outputs.

In this case, `multiext(config['genome'], ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")` generates the following six output files:

1. `config['genome'] + ".1.bt2"`
2. `config['genome'] + ".2.bt2"`
3. `config['genome'] + ".3.bt2"`
4. `config['genome'] + ".4.bt2"`
5. `config['genome'] + ".rev.1.bt2"`
6. `config['genome'] + ".rev.2.bt2"`

These files are the output of the "bowtie_index" rule, which builds a bowtie2 index from the reference genome.

*> conda: 'envs/bowtie2.yaml'*

means that the rule requires the use of a conda environment specified in the file "envs/bowtie2.yaml". The purpose of using a conda environment is to ensure that the necessary dependencies for the rule are installed and available. The yaml file specifies the packages and versions required for the environment. By using a conda environment, the pipeline can ensure that the correct dependencies are installed, even if they conflict with dependencies installed globally or used by other parts of the pipeline. This helps to prevent issues related to version conflicts and missing dependencies.

*> "threads: 2"*

specifies the number of threads to use in the processing step of the rule. In the rule "sam_to_bam", the line "threads: 2" means that the "samtools sort" command should use 2 threads in its processing. This allows the processing to be split across multiple cores, potentially reducing the processing time.

## Orginal provided workflow

```console
import pandas as pd
import re

configfile: 'config/config.yaml'

sample_df = (pd.read_csv(config['samples'],
    sep='\t',
    dtype={'sample_name':str, 'fastq':str})
    .set_index('sample_name', drop=False))



rule all:
    input:
        expand('results/01_called_variants/{sample}.vcf',
               sample=sample_df.sample_name)



rule bowtie_index:
    input:
        genome = config['genome']
    output:
        multiext(config['genome'],
                 ".1.bt2", ".2.bt2", ".3.bt2",
                 ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    conda: 'envs/bowtie2.yaml'
    shell:
        "bowtie2-build {input.genome} {input.genome} "



rule map_reads:
    input:
        idx = rules.bowtie_index.output,
        reads = lambda wildcards: sample_df.loc[wildcards.sample,
                                                 'fastq']
    output:
        temp('results/00_mapped_reads/{sample}.unsorted.sam')
    params:
        idx = config['genome']
    conda: 'envs/bowtie2.yaml'
    shell:
        'bowtie2 -x {params.idx} '
        '-U {input.reads} '
        '-S {output} '


rule sam_to_bam:
    input:
        rules.map_reads.output
    output:
        'results/00_mapped_reads/{sample}.bam'
    threads: 2
    conda: 'envs/htslib.yaml'
    shell:
        'samtools sort '
        '-@ {threads} '
        '-o {output} {input} '


rule index_bam:
    input:
        rules.sam_to_bam.output
    output:
        'results/00_mapped_reads/{sample}.bam.bai'
    conda: 'envs/htslib.yaml'
    shell:
        'samtools index {input} '


rule call_variants:
    input:
        rules.index_bam.output,
        aligned_reads = rules.sam_to_bam.output,
        genome = config['genome']
    output:
        'results/01_called_variants/{sample}.vcf'
    conda: 'envs/htslib.yaml'
    shell:
        'bcftools mpileup -Ou '
        '-f {input.genome} '
        '{input.aligned_reads} '
        '| bcftools call -m -v '
        '> {output} '

```

### Q8

* Get the above code up and running. Think about what the
  config/config.yaml file needs to contain. Note that you also need to set
  up the yaml files in the workflow/envs directory. Instructions for setting
  up environment files are available at https://snakemake.readthedocs.io/e
  n/stable/snakefiles/deployment.html#integrated-package-management.
  Let the htslib.yaml file contain version 1.11 of both samtools and
  bcftools. Note that both the conda-forge and bioconda channels are
  necessary to include in the environment file.

*config.yaml*

```console
samples: 'config/samples.tsv'
genome: 'resources/genome/yeastGenome.fa'
```

*samples.tsv*

```console
sample_name fastq
sample_1	resources/yeastReads.fastq
```

*bowtie2.yaml*

```console
channels:
- conda-forge
- bioconda

dependencies:
- bowtie2=2.4.2

```

*htslib.yaml*

```console
channels:
- conda-forge
- bioconda

dependencies:
- samtools=1.11
- bcftools=1.11
 
```

### Q9

* Use the workflow to rerun the variant calling. Use e.g. diff to ensure that
  the resulting vcf file looks similar to the one you generated in the variant
  calling exercise.

```console
$ snakemake -p --use-conda -j2
```

- Directory tree after running the workflow

```console
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
```

the bcftools headlines were different in the two vcf files

# Q10

I added newRead.fasq file into the same directory as the yeastReads.fastq file and modified the samples.tsv file to include the newReads.fastq file. the result directory tree is as follows:
the Q10VcSnakemakeArash.zip file contains all file necessary to run the snakemake workflow without intermediate and output files and also the .snakemake folder

```

.
├── config
│   ├── config.yaml
│   └── samples.tsv
├── image
│   └── README
├── README.md
├── resources
│   ├── genome
│   │   ├── yeastGenome.fa
│   │   ├── yeastGenome.fa.1.bt2
│   │   ├── yeastGenome.fa.2.bt2
│   │   ├── yeastGenome.fa.3.bt2
│   │   ├── yeastGenome.fa.4.bt2
│   │   ├── yeastGenome.fa.fai
│   │   ├── yeastGenome.fa.rev.1.bt2
│   │   └── yeastGenome.fa.rev.2.bt2
│   ├── newRead.fastq
│   └── yeastReads.fastq
├── results
│   ├── 00_mapped_reads
│   │   ├── newRead.bam
│   │   ├── newRead.bam.bai
│   │   ├── yeastReads.bam
│   │   └── yeastReads.bam.bai
│   └── 01_called_variants
│       ├── newRead.vcf
│       └── yeastReads.vcf
└── workflow
    ├── envs
    │   ├── bowtie2.yaml
    │   └── htslib.yaml
    └── Snakefile



```

## Q13

As optional part I chose Question13

* Write a new rule using the run or script property together with pandas
  or R to save the number of variants per scaffold to a tsv file. Beware that
  handling R dependencies in conda is not entirely straight forward.

After varinat calling pipline run command below in snakemake directory

```
snakemake  -p --use-conda -j2 --snakefile Snakefile2
```

which outputs tsv file/s in snakemake directory

> Snakemake2 rules

```
SAMPLE, = glob_wildcards('results/01_called_variants/{sample}.vcf')


rule all:
    input:
        expand('{sample}.tsv',sample=SAMPLE)

rule run_script:
    input:
       'results/01_called_variants/{sample}.vcf'
    output:
        '{sample}.tsv'
    shell:
        "python variants_per_scaffold.py {input} {output}"
```

> variants_per_scaffold.py script

```

import pandas as pd
import re
import sys

vcf_file = sys.argv[1]

# Read VCF into dataframe
df = pd.read_csv(vcf_file, sep='\t', comment='#', header=None)

# Extract scaffold name from the chromosome column (column 0)
df['scaffold'] = df[0].str.extract(r'(\w+)')

# Count number of variants per scaffold
scaffold_counts = df.groupby('scaffold').size().reset_index(name='count')

# Save scaffold counts to a TSV file
scaffold_counts.to_csv(sys.argv[2], sep='\t', index=False)
```

output tree for 2 fastq files:


```
.
├── config
│   ├── config.yaml
│   └── samples.tsv
├── image
│   └── README
├── README.md
├── resources
│   ├── genome
│   │   ├── yeastGenome.fa
│   │   ├── yeastGenome.fa.1.bt2
│   │   ├── yeastGenome.fa.2.bt2
│   │   ├── yeastGenome.fa.3.bt2
│   │   ├── yeastGenome.fa.4.bt2
│   │   ├── yeastGenome.fa.fai
│   │   ├── yeastGenome.fa.rev.1.bt2
│   │   └── yeastGenome.fa.rev.2.bt2
│   ├── newRead.fastq
│   └── yeastReads.fastq
├── results
│   ├── 00_mapped_reads
│   │   ├── newRead.bam
│   │   ├── newRead.bam.bai
│   │   ├── yeastReads.bam
│   │   └── yeastReads.bam.bai
│   └── 01_called_variants
│       ├── yeastReads2.vcf
│       └── yeastReads.vcf
├── Snakefile2
├── variants_per_scaffold.py
├── workflow
│   ├── envs
│   │   ├── bowtie2.yaml
│   │   └── htslib.yaml
│   └── Snakefile
├── yeastReads2.tsv
└── yeastReads.tsv

```
