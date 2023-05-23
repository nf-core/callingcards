# nf-core/callingcards: Output

## Introduction

The `callingcards` pipeline a workflow for both yeast and mammals callingcards
data. Both of these workflows executes similar steps. However, the details
and output are not identical. If you are interested in details of the step
or output, please ensure that you read the section corresponding to your
organism.

To run pre-configured tests on small data, suitable for a local computer,
choose your workflow type based on data type:

- Yeast

```bash
nextflow run nf-core/callingcards -profile test, <docker/singularity/institute>
```

- Mammals

```bash
nextflow run nf-core/callingcards -profile test_mammals, <docker/singularity/institute>
```

If you want to see the complete output, then set the following:

```
--save_genome_intermediate true --save_sequence_intermediate true --save_genome_intermediate true
```

### Output Overview

The output directory structure for both the Yeast and Mammals workflows is very similar. Notably, the

- [Preprocessing the Genome](#preprocessing-the-genome)

  - [Genome Masking](#genome-masking-bedtools)

  - [Adding additional sequences](#adding-additional-sequences-concatfasta)

  - [gtf format to bed format](#gtf-format-to-bed-format-gtf2bed)

  - [genome indexing](#aligner-indexing)

- [Sample Specific Output](#sample-specific-output)

  - [intermediate alignment results (selected aligner, samtools)](#alignment)

  - [Hops](#hops)

  - [QC modules (picard, rseqc, samtools)](#qc-modules-picard-rseqc-samtools)

  - [MultiQC](#multiqc)

- [Sequence](#sequence)

- [Pipeline info](#pipeline-info)

This pipeline provides two workflows, one to process yeast calling cards
data, and one to process mammal (mouse and human) data. The steps and output
are very similar, but there are subtle differences -- make sure that you are
reading the correct section.

### Preprocessing the Genome

If `save_genome_intermediate` is set to true (false by default), then the
`genome` subdirectory will be published in the `output_dir`.
Here is an example of the `genome` output. This was generated from the `test`
profile with `--save_genome_intermediate true`.

Note: the execution of any given run will be faster if the completely prepared
genome is passed through the appropriate parameters at runtime. If you don't
already have masked genome fasta file with the additional sequences added
(again, for yeast. For mammals this is not necessary), and the selected
aligner's index of this augmented genome, then run the pipeline once with
`--save_genome_intermediate true` and save the result for future use.

```bash
genome
├── bedtools
│   └── regions_mask.fa
├── bwamem2
│   └── bwamem2
│       ├── concat.fasta.0123
│       ...
├── concatfasta
│   └── concat.fasta
├── gtf2bed
│   └── genes.bed
└── samtools
    └── concat.fasta.fai
```

#### Genome Masking (bedtools)

A typical and recommended part of the yeast workflow, but not the mammal
workflow. When the user provides the `regions_mask` parameter (if you use
--genome R64-1-1 for yeast, this is set automatically by the pipeline),
the genome is hard masked by `bedtools maskfasta`.

#### Aligner indexing

For the chosen aligner, if that aligner's index for the genome.

#### Adding Additional Sequences (concatfasta)

A typical and recommended part of the yeast workfflow, but nto the mammal
workflow.

#### GTF format to Bed Format (gtf2bed)

This simply transforms the input GTF file to bed format.

#### Genome indexing (samtools)

This is the indexed genome -- indexing is performed after masking and,
critically, appending any additional sequences.

### Sample Specific Output

#### Alignment

After preprocessing the genome and raw reads, the processed reads are then
aligned to the processed genome. If `save_alignment_intermediate` is set to
false (false), then no `alignment` subdirectory will be created. Setting
`save_alignment_intermediate` to `true` would only be useful for debugging
purposes. An example of the output structure is:

```bash
alignment
├── bwamem2
│   ├── run_6177_T1_DAL80.bam
│   ├── ...
└── samtools
    ├── run_6177_T1_DAL80.bam.bai
    ├── ...
```

The alignment files which result from this step are input into the
[callingCardsTools](https://github.com/cmatKhan/callingCardsTools) counting
and sorting function, which determines if a read should be counted as a hop
and partitions the alignments into `passing` and `failing` bam files. It is
these partitioned files which are worth saving, not the intermediate
alignment files.

#### Hops

The yeast and mammals `hops` subdirectories contain slightly different QC
files, and contain slightly different intermediate subdirectories. It is
recommended that `save_alignment_intermediate` and `save_count_intermediate`
are both set to false (default) and not set to true except for debugging.
If you were to set both to true, you would see the following output:

**Yeast**

```bash
results/<sample>/hops
├── run_6177_T1_DAL80_failing_tagged.bam
├── run_6177_T1_DAL80_failing_tagged.bam.bai
├── run_6177_T1_DAL80_passing_tagged.bam
├── run_6177_T1_DAL80_passing_tagged.bam.bai
├── run_6177_T1_DAL80.qbed
├── run_6177_T1_DAL80_summary.tsv
├── ...
├── picard
│   ├── run_6177_T1_DAL80_failing_tagged.CollectMultipleMetrics.alignment_summary_metrics
│   ├── ...
├── rseqc
│   ├── run_6177_T1_DAL80_failing_tagged.read_distribution.txt
│   ├── ...
└── samtools
    ├── run_6177_T1_DAL80_failing_tagged.flagstat
    ├── run_6177_T1_DAL80_passing_tagged.flagstat
    ├── ...

```

- `*.passing/failing.bam(.bai)`

  - These are the alignment records, paritioned into those alignments which are
    counted as hops (passing) and those which are not (failing). the `.bai`
    file is the bam index

- `*.qbed`

  - A qbed file is a modified bed format file which describes Calling Cards
    transposon insertions.
    [See here for more information](https://cmatkhan.github.io/callingCardsTools/file_format_specs/qbed/).

- `*.summary.tsv`

  - This is the yeast alignment/hops level QC summary. Here is an example:

  ```raw
  status_decomp	count
  MAPQ	7
  MAPQ, FIVE_PRIME_CLIP	1
  NO_STATUS	17
  UNMAPPED	6
  ```

  where the first column is the quality status of a given alignment.
  `NO_STATUS` means passing, counted alignment. The second column are
  the number of alignment records in each category. Note that between any
  given sample, the levels of the `status_decomp` column may differ depending
  on how the reads classify.

**Mammals**

```bash
results/<sample>/hops
├── human_AY53-1_50_T1_passing_merged_sorted.bam
├── human_AY53-1_50_T1_passing_merged_sorted.bam.bai
├── human_AY53-1_50_T1_failing_merged_sorted.bam
├── human_AY53-1_50_T1_failing_merged_sorted.bam.bai
├── human_AY53-1_50_T1_aln_summary.tsv
├── human_AY53-1_50_T1_barcode_qc.tsv
├── human_AY53-1_50_T1.qbed
├── human_AY53-1_50_T1_srt_count.tsv
├── count
│   ├── ...
├── picard
│   ├── ...
└── samtools
    ├── ...

```

See the yeast hops section for an explanation of the `.bam(.bai)` and
`.qbed` files. Additionally, the `*_aln_summary.tsv` files are the same as
the yeast `summary.tsv` files -- see the yeast section above.

However, `srt_count.tsv` is unique to the mammals pipeline -- here is an
example of this file:

```raw
srt_type	count
single_srt	101
multi_srt	0
```

The first column,`srt_type` is either `single_srt` or `multi_srt`. The second
column `count`, describes how many, of the passing hops, of the insertions
have a single `srt` sequence, and how many have multiple `srt` sequences.

##### QC modules (picard, rseqc, samtools)

These are standard sequencing QC packages and the results are compiled by
MultiQC. See the [MultiQC](#multiqc) section for a discussion of how to
interpret these results for CallingCards experiments.

### Sequence

This subdirectory stores the result of the preprocessing of the raw reads. The
yeast and mammals workflows differ significantly in this case.

The **yeast** sequence directory looks like this:

```bash
sequence
├── run_6177_T1_null_r1_primer_summary.csv
├── run_6177_T1_null_r2_transposon_summary.csv
├── concatfastq
│   ├── run_6177_T1_DAL80_R1_concat.fastq
│   ├── run_6177_T1_DAL80_R2_concat.fastq
│   ├── ...
├── demultiplex
│   ├── barcode_qc_run_6177_T1_1.pickle
│   ├── DAL80_run_6177_T1_1_R1.fq
│   ├── DAL80_run_6177_T1_1_R2.fq
│   ├── ...
├── fastqcdemux
│   ├── run_6177_T1_DAL80_1_fastqc.html
│   ├── run_6177_T1_DAL80_1_fastqc.zip
│   ├── ...
├── fastqcraw
│   ├── ...
└── trimmomatic
    ├── ..
```

If `save_sequence_intermediate` is false (default, and recommended unless
debugging) then the output will be the summary files, and the fastqc files.

<details markdown="1">
<summary>Summary Files Description</summary>

- `*_primer_summary.csv`

  - This is a csv file which tallies the instances of the R2 transposon by
    R1 primer sequence (these two sequences together form the barcode which
    relates a transcription factor to a given read). This will look something
    like:

  ```raw
  tf,r1_primer_seq,r1_transposon_edit_dist,r2_transposon_edit_dist,restriction_ezyme,count
  MET31,TGATA,11,4,*,1
  INO2x1,CCTGC,13,7,*,1
  DAL80,CAACG,0,0,HinP1I,14
  DAL80,CAACG,0,0,*,2
  DAL80,CAACG,0,0,TaqAI,1
  DAL80,CAACG,0,0,Hpall,14
  DAL80,CAACG,0,7,TaqAI,1
  DAL80,CAACG,0,6,TaqAI,1
  DAL80,CAACG,0,6,HinP1I,1
  ```

- `*_transposon_summary.csv`
  - This is similar to the `_primer_summary.csv`, but in the reverse direction.
    This file instead tallies, for a given R2 transposon sequence, the number of
    `r1_primer_seq`s.

</details>

`demultiplex` stores the results of the partitioning of the reads by barcode
and `trimmomatic` stores the demultiplexed, trimmed reads. However, only
the `fastqc` directories will be saved if `save_sequence_intermediate` is
`false`, which is the default.

The **mammals** directory looks like this:

```bash
sequence/
├── fastqc
│   ├── human_AY53-1_50_T1_fastqc.html
│   └── human_AY53-1_50_T1_fastqc.zip
├── trimmomatic
│   ├── human_AY53-1_50_T1_1_barcoded_cropped.log
│   ├── human_AY53-1_50_T1_1_barcoded_cropped.SE.paired.trim.fastq.gz
│   ├── human_AY53-1_50_T1_1_barcoded_cropped.summary
│   ├── ...
└── umitools
    ├── human_AY53-1_50_T1_1_barcoded.umi_extract.fastq.gz
    ├── human_AY53-1_50_T1_1_barcoded.umi_extract.log
    ├── ...
```

where only `fastqc` will be present if `save_sequence_intermediate` is set to
`false` (default)

<details markdown="1">
<summary>FastQC Output Description</summary>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

<details markdown="1">
<summary>UMITools Extract Output Description</summary>

[UMItools extract](https://github.com/CGATOxford/UMI-tools) is used to extract barcode and other CallingCards specific sequences from
single or paired-end reads.

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

<details markdown="1">
<summary>Trimmomatic Output Description</summary>

[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) is used to trim reads after extracting known non-genomic
Calling Cards specific sequence

- `trimmomatic/`
  - `*.fastq.gz`: The fastq file, after the [`UMITools extract`](#umi-tools-extract)
  - `*.log`: Log file generated by tyrimmomatic

</details>


### MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

## Preprocessing Raw Reads

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

### Pipeline Info

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent
functionality for generating various reports relevant to the running and
execution of the pipeline. This will allow you to troubleshoot errors with the
running of the pipeline, and also provide you with other information such as
launch commands, run times and resource usage.

```bash
pipeline_info
├── execution_report_2023-05-21_13-43-36.html
├── execution_timeline_2023-05-21_13-43-36.html
├── execution_trace_2023-05-21_13-43-36.txt
├── pipeline_dag_2023-05-21_13-43-36.html
├── samplesheet.valid.csv
└── software_versions.yml
```

The `execution_report_...` is of particular interest if you wish to tune the
resource requests.

