//
// Align reads to a reference genome
// note that this can be parameterized -- could put $param.aligner
// in the include ... from ... path below
//

include { SAMTOOLS_SORT } from "${projectDir}/modules/nf-core/samtools/sort/main"
include { SAMTOOLS_INDEX } from "${projectDir}/modules/nf-core/samtools/index/main"
include { BWAMEM2_MEM    } from "${projectDir}/modules/nf-core/bwamem2/mem/main"
include { BWA_ALN        } from "${projectDir}/modules/nf-core/bwa/aln/main"
include { BWA_ALN_TO_BAM } from "${projectDir}/modules/local/bwa_aln_to_bam/main"
include { BOWTIE2_ALIGN    } from "${projectDir}/modules/nf-core/bowtie2/align/main"
include { BOWTIE_ALIGN    } from "${projectDir}/modules/nf-core/bowtie/align/main"

workflow ALIGN {
    take:
    reads           // channel: [ val(meta), [ reads ] ]
    bwamem2_index         // channel: file(fasta)
    bwa_aln_index
    bowtie2_index
    bowtie_index

    main:

    // NOTE: when adding an aligner, make sure the output is sorted and indexed.

    ch_versions = Channel.empty()
    ch_bam      = Channel.empty()

    if(params.aligner == 'bwamem2') {
        sort_bam = true
        BWAMEM2_MEM (
            reads,
            bwamem2_index,
            sort_bam
        )
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)

        SAMTOOLS_INDEX(
            BWAMEM2_MEM.out.bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

        // TODO figure out how to mix without using a tmp ch
        BWAMEM2_MEM.out.bam
            .join(SAMTOOLS_INDEX.out.bai)
            .set{ ch_tmp }

        ch_bam = ch_bam.mix(ch_tmp)

    } else if (params.aligner == 'bwaaln') {

        BWA_ALN (
            reads,
            bwa_aln_index
        )
        ch_versions = ch_versions.mix(BWA_ALN.out.versions)

        BWA_ALN_TO_BAM (
            reads.join(BWA_ALN.out.sai),
            bwa_aln_index
        )
        ch_versions = ch_versions.mix(BWA_ALN.out.versions)

        SAMTOOLS_INDEX(
            BWA_ALN_TO_BAM.out.bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

        // TODO figure out how to mix without using a tmp ch
        BWA_ALN_TO_BAM.out.bam
            .join(SAMTOOLS_INDEX.out.bai)
            .set{ ch_tmp }

        ch_bam = ch_bam.mix(ch_tmp)

    } else if (params.aligner == 'bowtie2') {
        sort_bam = true
        save_unaligned = false
        BOWTIE2_ALIGN (
            reads,
            bowtie2_index,
            save_unaligned,
            sort_bam
        )
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)

        SAMTOOLS_INDEX(
            BOWTIE2_ALIGN.out.bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

        // TODO figure out how to mix without using a tmp ch
        BOWTIE2_ALIGN.out.bam
            .join(SAMTOOLS_INDEX.out.bai)
            .set{ ch_tmp }

        ch_bam = ch_bam.mix(ch_tmp)

    } else if (params.aligner == 'bowtie') {

        BOWTIE_ALIGN (
            reads,
            bowtie_index
        )
        ch_versions = ch_versions.mix(BOWTIE_ALIGN.out.versions)

        SAMTOOLS_SORT(
            BOWTIE_ALIGN.out.bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

        SAMTOOLS_INDEX(
            SAMTOOLS_SORT.out.bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

        // TODO figure out how to mix without using a tmp ch
        SAMTOOLS_SORT.out.bam
            .join(SAMTOOLS_INDEX.out.bai)
            .set{ ch_tmp }

        ch_bam = ch_bam.mix(ch_tmp)

    } else {
        exit 1, "No aligner specified in params OR aligner: ${params.aligner} is not recognized. "
    }

    emit:
    bam       = ch_bam      // channel: [ val(meta), path(bam), path(bai) ]
    versions  = ch_versions // channel: [ versions.yml ]
}
