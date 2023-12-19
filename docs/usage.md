# nf-core/callingcards: Usage

<!-- ## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/callingcards/usage](https://nf-co.re/callingcards/usage) -->

<!-- > _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._ -->

## Introduction

Calling Cards experiments may be performed in both yeast and mammalian cells.
The appropriate workflow is selected with the `datatype` parameter. Suggested
default parameters for yeast and mammalian processing runs are provided through
the profiles [yeast](../conf/default_yeast.config) and
[mammal](../conf/default_mammal.config). These may be used by simply including
them with the `-profile` flag. See [Running the pipeline](#running-the-pipeline)
for examples of submission commands. See [`-profile`](#profile) for more
details on available profiles.

## Samplesheet input

You will need to create a samplesheet with information about the samples you
would like to analyse before running the pipeline. Use the `input` parameter to
specify its location. It has to be a comma-separated file with 4 columns, and
a header row as shown in the examples below.

**Note:** Currently, the mammals workflow supports only `fastq_1`. `fastq_2`
should simply be left blank. The yeast workflow requires both `fastq_1` and
`fastq_2`.

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using
the information provided in the samplesheet. The samplesheet can have as many
columns as you desire, however, there is a strict requirement for the first 4
columns to match those defined in the table below.

A final samplesheet file consisting of single end mammalian reads would look like so:

```console
sample,fastq_1,fastq_2,barcode_details
mouse_AY60-6_50,mouse/test_data/AY60-6_50k_downsampled_mouse.fastq.gz,,mouse/barcode_details.json
mouse_AY60-6_100,mouse/test_data/AY60-6_100k_downsampled_mouse.fastq.gz,,mouse/barcode_details.json
mouse_AY09-1_50_lowQuality,mouse/test_data/AY09-1_50k_downsampled_mouse_lowQuality.fastq.gz,,mouse/barcode_details.json
mouse_AY09-1_100_lowQuality,mouse/test_data/AY09-1_100k_downsampled_mouse_lowQuality.fastq.gz,,mouse/barcode_details.json

```

| Column            | Description                                                                                                                                                                            |
| ----------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`          | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1`         | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2`         | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `barcode_details` | Full path to the barcode details json file for a given sample.                                                                                                                         |

An [example samplesheet](../assets/human/input_samplesheet.csv) has been provided with the pipeline.

### Barcode Details

The barcode details json stores data which allows the pipeline to relate
sequence barcodes in the calling cards reads to a given transcription factor.

The file specificiations for both the yeast and mammals barcode details file
may be found [here](https://cmatkhan.github.io/callingCardsTools/file_format_specs/barcode_details/)

## Running the pipeline

The typical command for running the mammals workflow is as follows:

```bash
nextflow run nf-core/callingcards \
    -profile default_mammals,singularity \
    --input /path/to/your_samplesheet.csv \
    --fasta /path/to/genome.fa \
    --gtf /path/to/genes.gtf \
    --outdir results
```

A typical command for running the yeast workflow is as follows:

```bash
nextflow run nf-core/callingcards \
    -profile default_yeast,singularity \
    --input /path/to/your_samplesheet.csv \
    --outdir results
```

This will launch the pipeline with the specified profile(s). Note that the
pipeline will create the following files in your working directory:

```console
work            # (name configurable) Directory containing the nextflow working files
results         # (name configurable) Finished results
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

## General parameters

The following describes a selected set of parameters that are common to both
the yeast and mammalian workflows. For a full list of parameters, please see
the
[parameters section of the nf-core/callingcards site](https://nf-co.re/callingcards/dev/parameters)

- The `datatype` parameter accepts either `yeast` or `mammals` and determines
  which workflow to run.

- The `aligner` parameter accepts either `bwa`, `bwamem2`, `bowtie`, or `bowtie2`

- `split_fastq_by_size` or `split_fastq_by_part` controls how the fastq files
  are split for parallel processing. Set one or the other, not both.

- `min_mapq` sets the minimal mapping quality for reads to be considered
  'passing' in the hops counting stage. By default, this is `10`.

- `r1_crop` determines how much of R1 will be passed onto alignment. The read
  is cropped after extracting the non-genomic sequence.

- `save_genome_intermediate`, `save_sequence_intermediate`,
  and `save_alignment_intermediate` may be set to save intermediate files from
  each of the corresponding steps of the workflows.

## Mammals specific parameters

The mammals workflow requires that the following parameters be set. Note that
these parameters are set in the [default_mammals](../conf/default_mammals.config)
profile. But, you should confirm that these are correct for your data.

- `r1_bc_pattern` describes the barcode pattern that will be extracted by
  UMITools on R1.

## Yeast specific parameters

There are no yeast specific parameters -- rather the yeast specific steps are
entirely contained within the data provided in the barcode_details.json and
handled by [callingCardsTools](https://github.com/cmatKhan/callingCardsTools).
You should examine the [default_yeast](../conf/default_yeast.config)
configuration settings prior to using this profile.

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code
from GitHub and stores it as a cached version. When running the pipeline after
this, it will always use the cached version if available - even if the pipeline
has been updated since. To make sure that you're running the latest version of
the pipeline, make sure that you regularly update the cached version of the
pipeline:

```bash
nextflow pull nf-core/callingcards
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on
your data. This ensures that a specific version of the pipeline code and
software are used when you run your pipeline. If you keep using the same tag,
you'll be running the same version of the pipeline, even if there have been
changes to the code since.

First, go to the [nf-core/callingcards releases
page](https://github.com/nf-core/callingcards/releases) and find the latest
pipeline version - numeric only (eg. `1.3.1`). Then specify this when running
the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch
to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that
you'll know what you used when you look back in the future. For example, at the
bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter
files](#running-the-pipeline) to repeat pipeline runs with the same settings
without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for
academic publications), make sure to NOT include cluster specific paths to
files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline
parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give
configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the
pipeline to use software packaged using different methods (Docker, Singularity,
Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full
pipeline reproducibility, however when this is not possible, Conda is also
supported.
:::

The pipeline also dynamically loads configurations from
[https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it
runs, making multiple config profiles for various institutional clusters
available at run time. For more information and to see if your system is
available in these configs please see the [nf-core/configs
documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` -
the order of arguments is important! They are loaded in sequence, so later
profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all
software to be installed and available on the `PATH`. This is _not_ recommended,
since it can lead to different results on different machines dependent on the
computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/).
    Please only use Conda as a last resort i.e. when it's not possible to run
    the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud,
    or Apptainer.

#### Calling cards specific profiles

- `default_mammals`
  - A profile with suggested configuration for human and mouse data
- `default_yeast`
  - A profile with suggested configuration for human and mouse data
- `test`
  - A minimal test profile for the yeast workflow
- `test_mammals`
  - A minimal test profile for the mammalian workflow
- `test_full`
  - A minimal test profile for the full workflow -- mammals data

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from
any pipeline steps where the inputs are the same, continuing from where it got
to previously. For input to be considered the same, not only the names must be
identical but the files' contents as well. For more info about this parameter,
see [this blog
post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`.
Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command).
See the [nf-core website documentation](https://nf-co.re/usage/configuration)
for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for
most people and with most input data, you may find that you want to customise
the compute resources that the pipeline requests. Each step in the pipeline has
a default set of requirements for number of CPUs, memory and time. For most of
the steps in the pipeline, if the job exits with any of the error codes
specified
[here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18)
it will automatically be resubmitted with higher requests (2 x original, then 3
x original). If it still fails after the third attempt then the pipeline
execution is stopped.

To change the resource requests, please see the [max
resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning
workflow
resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources)
section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step
of the pipeline uses for a particular tool. By default nf-core pipelines use
containers and software from the [biocontainers](https://biocontainers.pro/) or
[bioconda](https://bioconda.github.io/) projects. However in some cases the
pipeline specified version maybe out of date.

To use a different container from the default container or conda environment
specified in a pipeline, please see the [updating tool
versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions)
section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a
particular tool used in pipeline. Fortunately, nf-core pipelines provide some
freedom to users to insert additional parameters that the pipeline does not
include by default.

To learn how to provide additional arguments to a particular tool of the
pipeline, please see the [customising tool
arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments)
section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if
you and others within your organisation are likely to be running nf-core
pipelines regularly and need to use the same settings regularly it may be a good
idea to request that your custom config file is uploaded to the
`nf-core/configs` git repository. Before you do this please can you test that
the config file works with your pipeline of choice using the `-c` parameter. You
can then create a pull request to the `nf-core/configs` repository with the
addition of your config file, associated documentation file (see examples in
[`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)),
and amending
[`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config)
to include your custom profile.

See the main [Nextflow
documentation](https://www.nextflow.io/docs/latest/config.html) for more
information about creating your own configuration files.

If you have any questions or issues please send us a message on
[Slack](https://nf-co.re/join/slack) on the [`#configs`
channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`. We recommend providing a compute `params.vm_type` of
`Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload
during the analysis. For a thorough list, please refer the [Azure Sizes for
virtual machines in
Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow
process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your
terminal so that the workflow does not stop if you log out of your session. The
logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a
detached session which you can log back into at a later time. Some HPC setups
also allow you to run nextflow within a cluster job submitted your job scheduler
(from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large
amount of memory. We recommend adding the following line to your environment to
limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
