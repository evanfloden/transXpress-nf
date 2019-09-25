params.reads = "$baseDir/data/reads/*_{1,2}.fq"
params.trimmomatic_adapters="$baseDir/data/adapters.fa"
params.outdir = "results"


log.info """\
 T R A N S X P R E S S - N F   P I P E L I N E
 ===================================
 reads        : ${params.reads}
 outdir       : ${params.outdir}
 """

Channel
    .fromFilePairs( params.reads, checkExists:true )
    .into{ read_pairs_ch; read_pairs_ch2 }


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


process trimmomatic {
    tag {"$reads[0]"+" and " +"$reads[1]"}

    input:
      set val(sample_id), file(reads) from read_pairs_ch2
      file 'adapters.fasta' from file(params.trimmomatic_adapters)

    output:
      set file("${reads[0]}.R1-P.qtrim.fastq.gz"), \
          file("${reads[1]}.R2-P.qtrim.fastq.gz") \
          into filteredPairedReads

    script:
      """
      trimmomatic PE -threads ${task.cpus} \
                     ${reads[0]} ${reads[1]} \
                     ${reads[0]}.R1-P.qtrim.fastq.gz \
                     ${reads[0]}.R1-U.qtrim.fastq.gz \
                     ${reads[1]}.R2-P.qtrim.fastq.gz \
                     ${reads[1]}.R2-U.qtrim.fastq.gz \
                     ILLUMINACLIP:adapters.fasta:2:30:10 \
                     SLIDINGWINDOW:4:5 \
                     LEADING:5 \
                     TRAILING:5 \
                     MINLEN:25 
      """
}
