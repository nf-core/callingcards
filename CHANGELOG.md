# nf-core/callingcards: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.0.0 - 2024-05-24

### Credits

Special thanks to the following for their reviews and assistance:

- [Maxime Garcia](https://github.com/maxulysse)
- [Friederike Hanssen](https://github.com/FriederikeHanssen)
- [Alyssa Briggs](https://github.com/alyssa-ab)
- [Sateesh Peri](https://github.com/sateeshperi)
- [James Fellows Yates](https://github.com/jfy133)

### `Fixed`

- workflow `datatype` option for mammals workflow corrected to `mammals`

- the modules running the callingCardsTools package were updated with the bioconda containers/conda paths. `callingCardsTools` updated to 1.2.0

- template update to v2.14.1 nf-core modules updated

- Nextflow `splitFastq` operator replaced by seqkit/split2 module

- multiqc output is reverted back to include all samples for the entire run, rather than included in each multiplexed libraries' subdirectories

### `Removed`

- local bwa_aln_to_bam local module in favor of nf-core fastq_align_bwa
  subworkflow

### `Added`

- nf-core/subworkflow fastq_align_bwa subworkflow

### `Changed`

- input fasta channel now creates a meta map w/ key id. This id is adjusted
  in the maskfasta and concatfasta module

- modules all updated to nf-core/modules current state

## 0.0.0dev - 2023-01-01

Initial release of nf-core/callingcards, created with the [nf-core](https://nf-co.re/) template.
