# nf-core/callingcards: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - 2023-10-31

### `Fixed`

- workflow `datatype` option for mammasl workflow corrected to `mammals`

- the modules running the callingCardsTools package were updated with the bioconda containers/conda paths. `callingCardsTools` updated to 1.2.0

- template update to v2.10. nf-core modules updated

- Nextflow `splitFastq` operator replaced by seqkit/split2 module

- multiqc output is reverted back to include all samples for the entire run, rather than included in each multiplexed libraries' subdirectories

## v0.0.0dev - 2023-01-01

Initial beta version
