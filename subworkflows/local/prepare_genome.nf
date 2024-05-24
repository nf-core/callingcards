//
// Index the reference genome and align reads
//

include { BEDTOOLS_MASKFASTA } from  "../../modules/nf-core/bedtools/maskfasta/main"
include { BWAMEM2_INDEX      } from  "../../modules/nf-core/bwamem2/index/main"
include { BWA_INDEX          } from  "../../modules/nf-core/bwa/index/main"
include { BOWTIE2_BUILD      } from  "../../modules/nf-core/bowtie2/build/main"
include { BOWTIE_BUILD       } from  "../../modules/nf-core/bowtie/build/main"
include { SAMTOOLS_FAIDX     } from  "../../modules/nf-core/samtools/faidx/main"
include { CONCATFASTA        } from  "../../modules/local/concatFasta/main.nf"
include { GTF2BED            } from  "../../modules/local/gtf2bed/main"

workflow PREPARE_GENOME {
    take:
    fasta         // [[id: 'identifier'], path(genome)]
    regions_mask  // [[id: 'identifier'], path(genome)]
    additional_fasta
    gtf

    main:

    ch_versions = Channel.empty()
    ch_bwamem2_index = Channel.empty()
    ch_bwa_index = Channel.empty()
    ch_bowtie2_index = Channel.empty()
    ch_bowtie_index = Channel.of("")

    // if a mask is provided, mask the fasta
    if(params.regions_mask){
        BEDTOOLS_MASKFASTA(
            regions_mask,
            fasta.map{meta, fasta -> fasta}
        )
        ch_versions = ch_versions.mix(BEDTOOLS_MASKFASTA.out.versions)
        ch_fasta = BEDTOOLS_MASKFASTA.out.fasta
            .map{meta, fasta -> [[id: 'masked_' + meta.id], fasta]}
    // else, just use the input fasta
    } else{
        ch_fasta = fasta
    }

    // if additional sequences are provided, append those to the fasta
    if (params.additional_fasta){
        CONCATFASTA(
            ch_fasta,
            additional_fasta
        )
        ch_versions = ch_versions.mix(CONCATFASTA.out.versions)

        ch_fasta = CONCATFASTA.out.fasta
            .map{meta, fasta -> [[id: 'concat_' + meta.id], fasta]}
    }

    GTF2BED(
        gtf
    )
    ch_versions = ch_versions.mix(GTF2BED.out.versions)

    if(params.aligner == 'bwamem2'){

        if (!params.bwamem2_index){
            BWAMEM2_INDEX ( ch_fasta )
            ch_bwamem2_index = BWAMEM2_INDEX.out.index
            ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
        } else{
            ch_bwamem2_index = Channel.of([[id: 'input_genome_index'],params.bwamem2_index]).collect()
        }

    } else if (params.aligner == 'bwa'){

        if (!params.bwa_index){
            BWA_INDEX ( ch_fasta )
            ch_bwa_index = BWA_INDEX.out.index
            ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
        } else{
            ch_bwa_index = Channel.of([[id: 'input_genome_index'],params.bwa_index]).collect()
        }

    } else if (params.aligner == 'bowtie2'){

        if (!params.bowtie2_index){
            BOWTIE2_BUILD ( ch_fasta )
            ch_bowtie2_index = BOWTIE2_BUILD.out.index
            ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
        } else{
            ch_bowtie2_index = Channel.of([[id: 'input_genome_index'],params.bowtie2_index]).collect()
        }

    } else if (params.aligner == 'bowtie'){

        if (!params.bowtie_index){
            BOWTIE_BUILD ( ch_fasta.map{meta, fasta -> fasta}  )
            ch_bowtie_index = BOWTIE_BUILD.out.index
            ch_versions = ch_versions.concat(BOWTIE_BUILD.out.versions)
        } else{
            ch_bowtie_index = Channel.of([params.bowtie_index]).collect()
        }

    } else {
        exit 1, "No aligner specified in params OR aligner: ${params.aligner} is not recognized. "
    }


    if (!params.fasta_index){
        SAMTOOLS_FAIDX (
            ch_fasta,
            [['id':null], []]
        )
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        ch_fasta_index = SAMTOOLS_FAIDX.out.fai
    } else {
        ch_fasta_index = Channel.of(["",params.fasta_index]).collect()
    }

    emit:
    fasta         = ch_fasta                // [val(meta), path(genome.fasta)]
    fai           = ch_fasta_index          // [val(meta), path(genome.fasta.fai)]
    genome_bed    = GTF2BED.out.bed         // [path(genome .bed)]
    bwamem2_index = ch_bwamem2_index        // [val(meta), index_data ]
    bwa_index     = ch_bwa_index            // [val(meta), index_data ]
    bowtie2_index = ch_bowtie2_index        // [val(meta), index_data ]
    bowtie_index  = ch_bowtie_index         // [index_data ]
    versions      = ch_versions             // channel: [ versions.yml ]
}
