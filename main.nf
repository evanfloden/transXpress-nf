params.reads = "$baseDir/data/reads/*_{1,2}.fq"
params.outdir = "results"

log.info """\
 T R A N S X P R E S S - N F   P I P E L I N E
 ===================================
 reads        : ${params.reads}
 outdir       : ${params.outdir}
 """

Channel
    .fromFilePairs( params.reads, checkExists:true )
    .set { read_pairs_ch }


process fastqc {
    tag "FASTQC on $sample_id"
    publishDir params.outdir

    input:
    set sample_id, file(reads) from read_pairs_ch

    output:
    file("fastqc_${sample_id}_logs") into fastqc_ch


    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}
