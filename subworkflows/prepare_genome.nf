//
// Index the reference genome and align reads
//

include { BWAMEM2_INDEX  } from  "${projectDir}/modules/nf-core/bwamem2/index/main"
include { BWA_INDEX      } from  "${projectDir}/modules/nf-core/bwa/index/main"
include { BOWTIE2_BUILD  } from  "${projectDir}/modules/nf-core/bowtie2/build/main"
include { BOWTIE_BUILD  } from  "${projectDir}/modules/nf-core/bowtie/build/main"
include { SAMTOOLS_FAIDX } from "${projectDir}/modules/nf-core/samtools/faidx/main"

workflow PREPARE_GENOME {
    take:
    fasta   // path(genome)

    main:

    ch_versions = Channel.empty()
    ch_bwamem2_index = Channel.empty()
    ch_bwa_aln_index = Channel.empty()
    ch_bowtie2_index = Channel.empty()
    ch_bowtie_index = Channel.of("")

    if(params.aligner == 'bwamem2'){

        if (!params.bwamem2_index){
            BWAMEM2_INDEX ( fasta.map{it -> ["", it]}  )
            ch_bwamem2_index = BWAMEM2_INDEX.out.index
            ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
        } else{
            ch_bwamem2_index = Channel.of(["",params.bwamem2_index]).collect()
        }

    } else if (params.aligner == 'bwaaln'){

        if (!params.bwaaln_index){
            BWA_INDEX ( fasta.map{it -> ["", it]}  )
            ch_bwa_aln_index = BWA_INDEX.out.index
            ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
        } else{
            ch_bwa_aln_index = Channel.of(["",params.bwa_aln_index]).collect()
        }

    } else if (params.aligner == 'bowtie2'){

        if (!params.bowtie2_index){
            BOWTIE2_BUILD ( fasta.map{it -> ["", it]}  )
            ch_bowtie2_index = BOWTIE2_BUILD.out.index
            ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
        } else{
            ch_bowtie2_index = Channel.of(["",params.bowtie2_index]).collect()
        }

    } else if (params.aligner == 'bowtie'){

        if (!params.bowtie_index){
            BOWTIE_BUILD ( fasta  )
            ch_bowtie_index = BOWTIE_BUILD.out.index
            ch_versions = ch_versions.concat(BOWTIE_BUILD.out.versions)
        } else{
            ch_bowtie_index = Channel.of([params.bowtie_index]).collect()
        }

    } else {
        exit 1, "No aligner specified in params OR aligner: ${params.aligner} is not recognized. "
    }


    if (!params.fasta_index){
        SAMTOOLS_FAIDX ( fasta.map{it -> ["", it]} )
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        ch_fasta_index = SAMTOOLS_FAIDX.out.fai
    } else {
        ch_fasta_index = Channel.of(["",params.fasta_index]).collect()
    }

    emit:
    fai           = ch_fasta_index          // [val(meta), path(*.fai)]
    bwamem2_index = ch_bwamem2_index        // [val(meta), index_data ]
    bwa_aln_index = ch_bwa_aln_index        // [val(meta), index_data ]
    bowtie2_index = ch_bowtie2_index        // [val(meta), index_data ]
    bowtie_index  = ch_bowtie_index         // [index_data ]
    versions      = ch_versions             // channel: [ versions.yml ]
}
