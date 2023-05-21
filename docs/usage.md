# nf-core/callingcards: Usage

<!-- ## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/callingcards/usage](https://nf-co.re/callingcards/usage) -->
<!-- ## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/callingcards/usage](https://nf-co.re/callingcards/usage) -->

<!-- > _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._ -->
<!-- > _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._ -->

## Introduction

Calling Cards experiments may be performed in both yeast and mammalian cells. The processing steps diverge, and this
divergence is controlled by setting certain parameters. Suggested default parameters for yeast and mammalian
processing runs are provided through the profiles [yeast](../conf/default_yeast.config) and [mammal](../conf/default_mammal.config). These may be used by simply including them with the `-profile` flag, for instance:

```
$ nextflow run callingcards/main.nf \
    -profile default_mammals,singularity \
    -c local.config \
    --input /path/to/samplesheet.csv
    --fasta /path/to/genome.fasta
    --output results_20220811
```

See [`-profile`](#profile) for more details on available profiles.

`local.config` is a [configuration](https://nf-co.re/usage/configuration) file which is the best place to put configuration settings such as what [executor](https://www.nextflow.io/docs/latest/executor.html)
you wish to use. If your institution already has a configuration profile on nf-core,
then you should use that profile in the `-profile` flag instead.
Calling Cards experiments may be performed in both yeast and mammalian cells. The processing steps diverge, and this
divergence is controlled by setting certain parameters. Suggested default parameters for yeast and mammalian
processing runs are provided through the profiles [yeast](../conf/default_yeast.config) and [mammal](../conf/default_mammal.config). These may be used by simply including them with the `-profile` flag, for instance:

```
$ nextflow run callingcards/main.nf \
    -profile default_mammals,singularity \
    -c local.config \
    --input /path/to/samplesheet.csv
    --fasta /path/to/genome.fasta
    --output results_20220811
```

See [`-profile`](#profile) for more details on available profiles.

`local.config` is a [configuration](https://nf-co.re/usage/configuration) file which is the best place to put configuration settings such as what [executor](https://www.nextflow.io/docs/latest/executor.html)
you wish to use. If your institution already has a configuration profile on nf-core,
then you should use that profile in the `-profile` flag instead.

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 4 columns, and a header row as shown in the examples below.
You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 4 columns, and a header row as shown in the examples below.

````console
```console
--input '[path to samplesheet file]'
````

Or, if you are using a file to save the run parameters (recommended), rather than submitting them on the command line,
then the file, with only the input set, would look like so:
Or, if you are using a file to save the run parameters (recommended), rather than submitting them on the command line,
then the file, with only the input set, would look like so:

````json
```json

{
  "input":"input_samplesheet.csv"
}
{
  "input":"input_samplesheet.csv"
}

````

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 4 columns to match those defined in the table below.
The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 4 columns to match those defined in the table below.

A final samplesheet file consisting of single end mammalian reads would look like so:
A final samplesheet file consisting of single end mammalian reads would look like so:

```console
sample,fastq_1,fastq_2,barcode_details
mouse_AY60-6_50,mouse/test_data/AY60-6_50k_downsampled_mouse.fastq.gz,,mouse/barcode_details.json
mouse_AY60-6_100,mouse/test_data/AY60-6_100k_downsampled_mouse.fastq.gz,,mouse/barcode_details.json
mouse_AY09-1_50_lowQuality,mouse/test_data/AY09-1_50k_downsampled_mouse_lowQuality.fastq.gz,,mouse/barcode_details.json
mouse_AY09-1_100_lowQuality,mouse/test_data/AY09-1_100k_downsampled_mouse_lowQuality.fastq.gz,,mouse/barcode_details.json

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

The barcode details json stores data which allows the pipeline to relate sequence barcodes in the calling cards reads to a given transcription factor.

#### **[Mammals](../assets/human/barcode_details.json)**

Yeast requires three fields: `indicies`, `components` and `tf_map`. Not that most likely,
you will only ever need to change the `tf_map` component.

| Fields         | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| -------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | --- |
| `batch`        | Optional. Might be a run, eg run_1234                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| `tf`           | Optional. Either ID or symbol of the TF                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| `r1`           | keys correspond to barcode component names. Each component is a map with keys `trim` which is boolean. Set to true to trim off this portion of the read (NOTE: not used in mammals pipeline. Instead controlled by UMITools). `index` specifies where this component occurs in the read                                                                                                                                                                                                                                                                                                                                                    |
| `r2`           | same as `r1`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| `components`   | Each key corresponds to a component in `r1` or `r2` and must be preceeded by `r1` or `r2` as appropriate. Each value is a map where the key `map` is required and lists the expected sequence(s) for that component. optional additional keys are `match_allowance` to allow mismatches > 0, `bam_tag` which is used to add the sequence extracted from the read to the aligned read, `match_type` which controls how the sequences are matched (default is `edit_distance` and appropriate in most circumstances), and `require` is a boolean which, when set to false, allows any number of mismatches in the component without penalty. |
| `insert_seq`   | A list of sequences expected to be present on the 5' end of the aligned read                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               | []  |
| `max_mismatch` | the maximum number of mismatches allowed in a barcode. By default this is equal to the sum of the mismatches allowed on each component. Set this to less than that sum to fail barcodes which pass component checks, but have more than x number of mismatches overall                                                                                                                                                                                                                                                                                                                                                                     | []  |

````json

{
    "batch": "",
    "tf": "",
    "r1": {
        "pb": {"trim": true,
               "index": [0,3]},

        "lrt1": {"trim": true,
                    "index": [3,28]},
        "srt": {"trim": true,
                "index":[28,32]},

        "lrt2": {"trim": true,
                    "index": [32,38]}
    },
    "r2":{},
    "components": {

        "r1_pb":    {"map":["TAG"],
                     "match_allowance": 0,
                     "bam_tag": "PB"},

        "r1_lrt1":  {"map": ["CGTCAATTTTACGCAGACTATCTTT"],
                     "match_type": "edit_distance",
                     "match_allowance": 0,
                     "require": true,
                     "bam_tag": "L1"},
        "r1_srt":   {"map": ["CTAG", "CAAC", "CTGA", "GCAT", "GTAC", "CACA", "TGAC", "GTCA",
                             "CGAT", "CTCT", "GAAG", "TCGA", "CATG", "GTTG", "CTTC", "GCTA",
                             "GAGA", "GTGT", "CGTA", "TGGT", "GGAA", "ACAC", "TCAG", "TTGG",
                             "CAGT", "TTTT"],
                     "match_type": "edit_distance",
                     "match_allowance": 0,
                     "require": true,
                     "bam_tag": "ST"},
        "r1_lrt2":  {"map": ["GGTTAA"],
                     "match_type": "edit_distance",
                     "match_allowance": 0,
                     "require": true,
                     "bam_tag": "L2"}
    },
    "insert_seq": ["TTAA"],
    "max_mismatch": 0
}

| Column         | Description                                                                                                                                                                            |
|----------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `sample`          | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1`         | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2`         | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".
| `barcode_details` | Full path to the barcode details json file for a given sample. |

An [example samplesheet](../assets/human/input_samplesheet.csv) has been provided with the pipeline.

### Barcode Details

The barcode details json stores data which allows the pipeline to relate sequence barcodes in the calling cards reads to a given transcription factor.

#### __[Mammals](../assets/human/barcode_details.json)__

Yeast requires three fields: `indicies`, `components` and `tf_map`. Not that most likely,
you will only ever need to change the `tf_map` component.

| Fields         | Description                                                                                                                                                                            |
|----------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `batch`        | Optional. Might be a run, eg run_1234 |
| `tf`           | Optional. Either ID or symbol of the TF |
| `r1`           | keys correspond to barcode component names. Each component is a map with keys `trim` which is boolean. Set to true to trim off this portion of the read (NOTE: not used in mammals pipeline. Instead controlled by UMITools). `index` specifies where this component occurs in the read |
| `r2`           | same as `r1` |
| `components`   | Each key corresponds to a component in `r1` or `r2` and must be preceeded by `r1` or `r2` as appropriate. Each value is a map where the key `map` is required and lists the expected sequence(s) for that component. optional additional keys are `match_allowance` to allow mismatches > 0, `bam_tag` which is used to add the sequence extracted from the read to the aligned read, `match_type` which controls how the sequences are matched (default is `edit_distance` and appropriate in most circumstances), and `require` is a boolean which, when set to false, allows any number of mismatches in the component without penalty.|
| `insert_seq`   | A list of sequences expected to be present on the 5' end of the aligned read |[]
| `max_mismatch` | the maximum number of mismatches allowed in a barcode. By default this is equal to the sum of the mismatches allowed on each component. Set this to less than that sum to fail barcodes which pass component checks, but have more than x number of mismatches overall |[]

```json

{
    "batch": "",
    "tf": "",
    "r1": {
        "pb": {"trim": true,
               "index": [0,3]},

        "lrt1": {"trim": true,
                    "index": [3,28]},
        "srt": {"trim": true,
                "index":[28,32]},

        "lrt2": {"trim": true,
                    "index": [32,38]}
    },
    "r2":{},
    "components": {

        "r1_pb":    {"map":["TAG"],
                     "match_allowance": 0,
                     "bam_tag": "PB"},

        "r1_lrt1":  {"map": ["CGTCAATTTTACGCAGACTATCTTT"],
                     "match_type": "edit_distance",
                     "match_allowance": 0,
                     "require": true,
                     "bam_tag": "L1"},
        "r1_srt":   {"map": ["CTAG", "CAAC", "CTGA", "GCAT", "GTAC", "CACA", "TGAC", "GTCA",
                             "CGAT", "CTCT", "GAAG", "TCGA", "CATG", "GTTG", "CTTC", "GCTA",
                             "GAGA", "GTGT", "CGTA", "TGGT", "GGAA", "ACAC", "TCAG", "TTGG",
                             "CAGT", "TTTT"],
                     "match_type": "edit_distance",
                     "match_allowance": 0,
                     "require": true,
                     "bam_tag": "ST"},
        "r1_lrt2":  {"map": ["GGTTAA"],
                     "match_type": "edit_distance",
                     "match_allowance": 0,
                     "require": true,
                     "bam_tag": "L2"}
    },
    "insert_seq": ["TTAA"],
    "max_mismatch": 0
}


````

```

## Running the pipeline

The typical command for running the pipeline is as follows:

```

$ nextflow run callingcards/main.nf \
 -profile default_mammals,singularity \
 -c local.config \
 --input /path/to/your_samplesheet.csv
--fasta /path/to/genome.fa

````

__Note__: you can also save your input parameters in a json format. For
instance, in the case of the run command above, you could create a
file called `params.json` which looks like so:

```json
{
    "fasta":"/path/to/genome.fa",
    "input":"/path/to/your_samplesheet.csv"
}

````

in which case the `run` command would look like:

```
$ nextflow run callingcards/main.nf \
    -profile default_mammals,singularity \
    -c local.config \
    --input /path/to/your_samplesheet.csv
    --fasta /path/to/genome.fa
```

**Note**: you can also save your input parameters in a json format. For
instance, in the case of the run command above, you could create a
file called `params.json` which looks like so:

```json
{
  "fasta": "/path/to/genome.fa",
  "input": "/path/to/your_samplesheet.csv"
}
```

in which case the `run` command would look like:

```bash
$ nextflow run callingcards \
    -profile default_mammals,singularity \
    -c local.config \
    -params-file path/to/params.json

$ nextflow run callingcards \
    -profile default_mammals,singularity \
    -c local.config \
    -params-file path/to/params.json

```

This will launch the pipeline with the `default_mammal` pipeline settings. It will use the `singularity` profile, and
set user specific system configuration (eg executor settings) in a file called local.config. Parameters such as `input`,
the path to the samplesheet, `output`, and other run parameters will be set in the params.json file.
See below for more information about profiles.

This will launch the pipeline with the `default_mammal` pipeline settings. It will use the `singularity` profile, and
set user specific system configuration (eg executor settings) in a file called local.config. Parameters such as `input`,
the path to the samplesheet, `output`, and other run parameters will be set in the params.json file.
See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

````console
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
```console
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
````

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> ⚠️ Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
> The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/callingcards -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
input: 'data'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/callingcards
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/callingcards releases page](https://github.com/nf-core/callingcards/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.
If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

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
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
- `default_mammals`
  - A profile with suggested configuration for human and mouse data
- `test`
  - A profile with a complete configuration for human data for automated testing
  - Mouse data would run through the same steps, hence `default_mammals` is
    appropriate for both human and mouse
  - Includes links to test data -- no other parameters necessary

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

```console
[62/149eb0] NOTE: Process `RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

Command error:
    .command.sh: line 9:  30 Killed    STAR --genomeDir star --readFilesIn WT_REP1_trimmed.fq.gz --runThreadN 2 --outFileNamePrefix WT_REP1. <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

#### For beginners

A first step to bypass this error, you could try to increase the amount of CPUs, memory, and time for the whole pipeline. Therefor you can try to increase the resource for the parameters `--max_cpus`, `--max_memory`, and `--max_time`. Based on the error above, you have to increase the amount of memory. Therefore you can go to the [parameter documentation of rnaseq](https://nf-co.re/rnaseq/3.9/parameters) and scroll down to the `show hidden parameter` button to get the default value for `--max_memory`. In this case 128GB, you than can try to run your pipeline again with `--max_memory 200GB -resume` to skip all process, that were already calculated. If you can not increase the resource of the complete pipeline, you can try to adapt the resource for a single process as mentioned below.

#### Advanced option on process level

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process. The quickest way is to search for `process STAR_ALIGN` in the [nf-core/rnaseq Github repo](https://github.com/nf-core/rnaseq/search?q=process+STAR_ALIGN).
We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so, based on the search results, the file we want is `modules/nf-core/star/align/main.nf`.
If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_high`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L9).
The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements.
The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L33-L37) which in this case is defined as 72GB.
Providing you haven't set any other standard nf-core parameters to **cap** the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `STAR_ALIGN` process failure by creating a custom config file that sets at least 72GB of memory, in this case increased to 100GB.
The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: STAR_ALIGN {
        memory = 100.GB
    }
}
```

> **NB:** We specify just the process name i.e. `STAR_ALIGN` in the config file and not the full task name string that is printed to screen in the error message or on the terminal whilst the pipeline is running i.e. `RNASEQ:ALIGN_STAR:STAR_ALIGN`. You may get a warning suggesting that the process selector isn't recognised but you can ignore that if the process name has been specified correctly. This is something that needs to be fixed upstream in core Nextflow.

### Updating containers (advanced users)

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

   - For Docker:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Singularity:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Conda:

     ```nextflow
     process {
         withName: PANGOLIN {
             conda = 'bioconda::pangolin=3.0.5'
         }
     }
     ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

````console
```console
NXF_OPTS='-Xms1g -Xmx4g'
````
