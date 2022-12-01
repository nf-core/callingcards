# ![nf-core/callingcards-mammals](docs/images/nf-core-callingcards-mammals_logo_light.png#gh-light-mode-only) ![nf-core/callingcards-mammals](docs/images/nf-core-callingcards-mammals_logo_dark.png#gh-dark-mode-only)

[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/callingcards-mammals/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/callingcards-mammals)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23callingcards-mammals-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/callingcards-mammals)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->

**nf-core/callingcards-mammals** is a bioinformatics best-practice analysis pipeline for An automated processing pipeline for mammalian bulk calling cards experiments.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

The testing data is currently provided in the `assets` directory of the pipeline code repository. Instructions on running the tests area available in Usage.

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/callingcards/results).

## Pipeline summary

1. Prepare Reads
    1. Extract barcodes ([`UMItools`](https://github.com/CGATOxford/UMI-tools))
    1. Trim, and reduce to only R1 depending on user input ([`Trimmomatic`](http://www.usadellab.org/cms/?page=trimmomatic))
1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
1. Prepare the Genome
    1. Samtools faidx and aligner indicies
1. Alignment
    1. One of: [`bwamem2`](https://github.com/bwa-mem2/bwa-mem2),[`bwa`](https://bio-bwa.sourceforge.net/bwa.shtml),[`bowtie2`](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml),[`bowtie`](https://bowtie-bio.sourceforge.net/index.shtml)
1. Process Alignments
    1. Extract alignment QC metrics ([`Samtools`](https://www.htslib.org/), [`preseq`](http://smithlabresearch.org/software/preseq/))
    1. Add Read Group and tags to alignment files ([custom script](https://github.com/cmatKhan/pycallingcards/tree/raw_processing/pycallingcards/raw_processing))
    1. Quantify transposon hops ([custom script](https://github.com/cmatKhan/pycallingcards/tree/raw_processing/pycallingcards/raw_processing))
1. Present QC for raw read and alignment metrics ([`MultiQC`](http://multiqc.info/))

## Quick Start

Note that more detailed instructions are available in [usage](docs/usage.md).
Specific instructions for running this on [HTCF are provided here](docs/htcf_specific.md)

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

1. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

1. Git clone the repo and test it on a minimal data set designed to confirm that the pipeline will complete

   ```bash
   $ mkdir calling_card_output

   $ cd calling_card_output

   $ git clone https://github.com/cmatKhan/callingcards-mammals.git

   # you can replace singularity below with any one of the dependency managers
   # There are also other test profiles -- see usage for more detail
   $ nextflow run callingcards-mammals/main.nf -profile test_human,singularity
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`test_yeast` and `singularity` in the example command above). You can chain multiple config profiles in a comma-separated string, as demonstrated.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

1. Start running your own analysis!

   ```bash
   nextflow run callingcards-mammals/main.nf \
        -params-file params.json \
        -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> \
        # possibly more config settings for your environment
        -c local.config
   ```

The `params.json` file is described in [usage](docs/usage.md)
Configuration is discussed in [Pipeline configuration](https://nf-co.re/usage/configuration) and
in the [configuration section of the nextflow documentation](https://www.nextflow.io/docs/latest/config.html)

## Documentation

The nf-core/callingcards pipeline comes with documentation about the pipeline [usage](docs/usage.md), [parameters](docs/parameters.md) and [output](docs/output.md).

## Credits

nf-core/callingcards was written by Chase Mateusiak. It was adapted from scripts written by: __fill this in__

We thank the following people for their extensive assistance in the development of this pipeline: __fill this in__

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#callingcards` channel](https://nfcore.slack.com/channels/callingcards) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/callingcards for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
