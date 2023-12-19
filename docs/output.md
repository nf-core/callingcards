# nf-core/callingcards: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots
are taken from the MultiQC report, which summarises results at the end of the
pipeline.

In the output directory of a given run of this pipeline, there will be a
subdirectory for each sample in the [samplesheet](./usage.md#samplesheet-input)
where the directory name is the sample name (first field in the samplesheet).

For example, if you are running **mammals** data and your samplesheet has two
samples, like so:

```raw
sample,fastq_1,fastq_2,barcode_details
AAV9_cortex,data/AAV9_1-1_cortex_L004_R1.fastq.gz,,barcode_details.json
AAV9_hindbrain,data/AAV9_1-1_hindbrain_L004_R1.fastq.gz,,barcode_details.json
```

Then the output directory would have the following structure:

```raw
example_results/
├── AAV9_cortex
│   ├── hops
│   └── sequence
├── AAV9_hindbrain
│   ├── hops
│   └── sequence
├── multiqc
└── pipeline_info

```

If you are running **yeast** data, then the samplesheet would look the same,
except fastq_2 must be present. The output directory will be structured
similarly, except that within each sample subdirectory, there will be multiple
qbeds/qc files for each of the TFs included in the multiplexed library. Mammals
libraries, on the other hand, are not multiplexed and contain data for only 1
TF.

For either `datatype`, you can include the parameters
`--save_genome_intermediate`, `--save_sequence_intermediate` and/or
`--save_alignment_intermeidate` to save the intermediate files generated at
each processing stage.

The directories listed below will be created in the results directory after the
pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

1. Prepare the Genome
   - [Mask Fasta](#mask-genome)
   - [Concat Fasta](#concatenate-additional-sequences)
   - [GTF2BED](#gtf2bed)
   - [Aligner Index](#aligner-index)
   - [Genome Index](#genome-index)
2. Prepare the Reads
   - [FastQC](#fastqc)
   - [Divide Fastq](#divide-fastq)
   - [UMITools Extract](#umitools-extract) (mammals only)
   - [Demultiplex](#demultiplex) (yeast only)
   - [Concat Fastqs](#concat-fastq) (yeast only)
   - [Concat QC](#concat-qc) (for mammals, this occurs in the Align step)
   - [Trimmomatic](#trimmomatic)
3. Align
   - [Align](#align)
4. Count Hops
   - [Count Hops](#count-hops)
   - [Alignment QC and Metrics](#alignment-metrics)
5. [Present QC data](#multiqc)
6. [Pipeline information](#pipeline-information)

### Mask Genome

This step does not output, except in the `work-dir`.
However, if a bed file is provided to the `regions_mask` parameter, then it
will be used to mask the genome with `bedtools maskfasta`. The masked genome
will be used for the rest of the pipeline.

### Concatenate Additional Sequences

This step does not output, except in the `work-dir`.
However, if a fasta file is provided to the the `additional_sequence`
parameter, it will be appended to the (possibly masked) genome and used for all
subsequent steps.

### GTF2BED

Translation of the gtf file to a bed file.
This is only output if `save_genome_intermediate` is `true`

<details markdown="1">
<summary>Output files</summary>

- `genome/gtf2bed`
  - `<gtf_name>.bed`: GTF file in bed format

</details>

### Aligner Index

The index produced by the aligner of your choice. This will only be saved if
`save_genome_intermediate` is true

<details markdown="1">
<summary>Output files</summary>

- `genome/<aligner>`
  - `<aligner index output>`: The output of the aligner index command

</details>

### Genome Index

Samtools index of the fasta file (after appending any additional sequences)

<details markdown="1">
<summary>Output files</summary>

- `genome/samtools`
  - `<genome>.fa[asta].fai`: Genome index (includes additional sequences, if they exist).

</details>

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The yeast workflow runs FastQC both before and after processing (including demultiplexing) the reads. The
> mammalian workflow runs FastQC only on un-processed reads

### Divide Fastq

seqkit split2 is used to divide each of the samples' fastq files into parts
for parallel processing. For the yeast pipeline, these parts are concatenated
back together by the respective TFs while for mammals, the divided parts are
aligned and processed in parallel

### UMItools Extract

<details markdown="1">
<summary>Output files</summary>

- `umitools/`

  - `*.fastq.gz`: The input fastq files with the specified pattern extracted and placed in the
    header. For example, for single-end reads where the `r1_bc_pattern` is
    `NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN`, the result of this step would turn a FASTQ record like this:

  ```raw
    @A00118:503:HNGFYDSX3:4:1101:2844:1000 1:N:0:CCTACCAAAG+TATGTTTATG
    TAGCGTCAATTTTACGCAGACTATCTTTCTGAGGTTAAAGTCTACCCACTTATTGAATTGTATATGTGTGATTCAGTTTTGCCAATGATGCTTTCCAGTCTCTAAAACGGATTTATCCTGGGTAACAAGAGCCTTGCACATTTTGATGACT
    +
    FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFF

  ```

  to this:

  ```raw
    @A00118:503:HNGFYDSX3:4:1101:2844:1000_TAGCGTCAATTTTACGCAGACTATCTTTCTGAGGTTAA 1:N:0:CCTACCAAAG+TATGTTTATG
    AGTCTACCCACTTATTGAATTGTATATGTGTGATTCAGTTTTGCCAATGATGCTTTCCAGTCTCTAAAACGGATTTATCCTGGGTAACAAGAGCCTTGCACATTTTGATGACT
    +
    FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFF
  ```

  Where the barcode immediately follows the first space delimited value in the id line and is appended with an `_`

  - `*.log`: Log file generated by the UMI-tools `extract` command.

</details>

[UMItools extract](https://github.com/CGATOxford/UMI-tools) is used to extract barcode and other CallingCards specific sequences from
single or paired-end reads.

### Demultiplex

This step is only run for the yeast pipeline. It calls the
[callingCardsTools](https://github.com/cmatKhan/callingCardsTools) utility
`parse fastq` and uses the `barcod_details` json file to demultiplex the reads
according to the TF barcodes. This step will typically be performed on
split fastq files (for parallel processing) and not output unless
`save_sequence_itermediate` is set.

### Concat Fastq

This step is only run for the yeast pipeline. It concatenates the splits of
each demultiplexed TF fastq file so that there is a single fastq file for each
TF for each sample.

### Concat QC

This is performed in both the yeast and mammals workflows, but at different
points. For yeast, this is performed in the
[Prepare Reads](../subworkflows/local/yeast/prepare_reads.nf) subworkflow while
in mammals it is performed after alignment in the
[Process Alignment](../subworkflows/local/mammals/process_alignments.nf)
subworkflow.

<details markdown="1">
<summary>Output files</summary>

- `sequence/<sample>`
  - `fastqcdemux`: The trimmed fastq file
  - `fastqcraw`: Log file generated by tyrimmomatic
  - `*_r1_primer_summary`: Given a correct R1 barcode, this file describes how many reads go to each R2 barcode tallied by edit distance. So, the entry
    `MET31,TGATA,11,4,*,1` means that for the MET31 R1 sequence `TGATA` which has an R1 transposon sequence
    with an edit distance of 11 and R2 TF barcode sequence
    with an edit distance of 4 and and unknown restriction
    enzyme, there was 1 read.
  - `*_r2_transposon_summary.csv`: Similar to above, but
    rather than tallying by correct R1 sequences, this tallies by correct R2 TF barcode sequences.

</details>

### Trimmomatic

<details markdown="1">
<summary>Output files</summary>

- `trimmomatic/`
  - `*.fastq.gz`: The trimmed fastq file
  - `*.log`: Log file generated by tyrimmomatic

</details>

[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) is used to trim reads after extracting known non-genomic
Calling Cards specific sequence

### Align

The user may select any one of the following:

- [bwa aln](https://github.com/lh3/bwa)
- [bwa mem2](https://github.com/bwa-mem2/bwa-mem2)
- [bowtie](https://bowtie-bio.sourceforge.net/index.shtml)
- [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

If `save_intermediate` is set, and the selected Aligner generates an index,
the index will be saved at the main level of the `outdir` in a directory called
`genome`. In each sample subdirectory, there will also be a directory
named by the name of the selected aligner which will store the aligner's output.

<details markdown="1">
<summary>Output files</summary>

- `genome/`

  - `<ALIGNER>`: this will store the selected aligners index, if `save_intermediate` is set.

- `<ALIGNER>/`
  - `<SAMPLE>.sorted.bam`: Coordinate sorted alignment file generated by the user-selected aligner
  - `<SAMPLE>.sorted.bam.bai`: the BAI index for the alignment file
  </details>

### Count Hops

<details markdown="1">
<summary>Output files</summary>

- `hops/`
  - `*_passing.bam/.bai`: the subset of alignments which are considered passing, countable hops,
    and its BAI index
  - `*_failing.bam/.bai`: the subset of alignments which are considered uncountable
  - `*.qbed`: A qBed format file which quantifies the number of hops at a given
    coordinate in the genome
  - `*_summary.tsv`: A tally of the number of reads which passed and failed by failure status
  - `*_barcode_qc.tsv`: (mammals only) A tally of the number of reads by barcode components
  - `*_srt_count.tsv`: (mammals only) A tally of reads with a single and multi SRT sequence

</details>

### Alignment Metrics

Picard, RSeQC and Samtools provide alignment level quality control metrics.
Qualifying alignments are counted as 'hops' of a given transcription factor
and those hops are quantified in a qBed file. Hop level QC metrics are also
generated.

#### Picard

Picard CollectMultipleMetrics gathers multiple QC metrics from alignments files. These are performed on
the post-processed alignments files, after they have
been partitioned into passing and failing alignments.

<details markdown="1">
<summary>Output files</summary>

- `hops/picard/`
  - `*.alignment_summary_metrics`
  - `.base_distribution_by_cycle_metrics/.pdf`
  - `*.quality_by_cycle_metrics/.pdf`
  - `*quality_distribution_metrics/.pdf`
  - `*.read_length_histogram.pdf`

</details>

#### SAMtools

<details markdown="1">
<summary>Output files</summary>

- `alignment/`
  - `<SAMPLE>.sorted.bam.bai`: the BAI index for the alignment file
- `hops/samtools/`
  - SAMtools `<SAMPLE>.sorted.bam.flagstat`, `<SAMPLE>.sorted.bam.idxstats` and `<SAMPLE>.sorted.bam.stats` files generated from the alignment files.

</details>

The original BAM files generated by the selected alignment algorithm are further processed with [SAMtools](http://samtools.sourceforge.net/) to sort them by coordinate, for indexing, as well as to generate read mapping statistics.

#### RSeQC

[RSeQC](<(http://rseqc.sourceforge.net/)>) is a package of scripts designed to evaluate the quality of RNA-seq data. This pipeline runs several, but not all RSeQC scripts. You can tweak the supported scripts you would like to run by adjusting the `--rseqc_modules` parameter which by default will run all of the following: `bam_stat.py`, `inner_distance.py`, `infer_experiment.py`, `junction_annotation.py`, `junction_saturation.py`,`read_distribution.py` and `read_duplication.py`.

The majority of RSeQC scripts generate output files which can be plotted and summarised in the MultiQC report.

##### Infer experiment

<details markdown="1">
<summary>Output files</summary>

- `hops/rseqc`
  - `*.infer_experiment.txt`: File containing fraction of reads mapping to given strandedness configurations.

</details>

A number of RSeQC modules are available. Only the read distribution

RSeQC documentation: [infer_experiment.py](http://rseqc.sourceforge.net)

![MultiQC - Strand check table](images/mqc_strand_check.png)

##### Read distribution

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/rseqc/read_distribution/`
  - `*.read_distribution.txt`: File containing fraction of reads mapping to genome feature e.g. CDS exon, 5’UTR exon, 3’ UTR exon, Intron, Intergenic regions etc.

</details>

This tool calculates how mapped reads are distributed over genomic features. A good result for a standard RNA-seq experiments is generally to have as many exonic reads as possible (`CDS_Exons`). A large amount of intronic reads could be indicative of DNA contamination in your sample but may be expected for a total RNA preparation.

RSeQC documentation: [read_distribution.py](http://rseqc.sourceforge.net/#read-distribution-py)

##### BAM stat

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/rseqc/bam_stat/`
  - `*.bam_stat.txt`: Mapping statistics for the BAM file.

</details>

This script gives numerous statistics about the aligned BAM files.

MultiQC plots each of these statistics in a dot plot. Each sample in the
project is a dot - hover to see the sample highlighted across all fields.

RSeQC documentation: [bam_stat.py](http://rseqc.sourceforge.net/#bam-stat-py)

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
