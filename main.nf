params.reads = "$baseDir/data/reads/*_{1,2}.fq"
params.trimmomatic_adapters="$baseDir/data/adapters.fa"
params.outdir = "results"
params.rnaSPAdes_params = ""


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
      set val(sample_id),\
          file("${reads[0]}.R1-P.qtrim.fastq.gz"), \
          file("${reads[1]}.R2-P.qtrim.fastq.gz") \
          into filteredPairedReads_ch1, \
               filteredPairedReads_ch2

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

filteredPairedReads_ch1
  .collectFile { it[0]+'\t'+"${it[0]}-rep1"+'\t'+it[1].name+'\t'+it[2].name+'\n'}
  .set {rnaSPAdes_samples_ch}


process rnaSPAdes_yaml {

    input:
      file rnaSPAdes_sample from rnaSPAdes_samples_ch

    output:
      file "rnaSPAdes_samples.yaml" into rnaSPAdes_yaml_ch

script:
    """
    #!/usr/bin/env python
    import re
    import os
    import os.path
    import pprint
    sample_list = []
    with open("${rnaSPAdes_sample}", "r") as input_handle:
      for line in input_handle:
        row = re.split("[\t ]+", line)
        if (len(row) > 3): # paired reads
          paired_dict = {}
          paired_dict['orientation'] = 'fr'
          paired_dict['type'] = 'paired-end'
          f_reads = row[2].strip()
          r_reads = row[3].strip()
          paired_dict['left reads'] = [f_reads]
          paired_dict['right reads'] = [r_reads]
          sample_list.append(paired_dict)
        if (len(row) == 3): # unpaired reads
          unpaired_dict = {}
          unpaired_dict['type'] = 'single'
          u_reads = row[2].strip()
          assert os.path.isfile(u_reads)
          unpaired_dict['single reads'] = [u_reads]
          sample_list.append(unpaired_dict)
    with open("rnaSPAdes_samples.yaml", "w") as output_handle:
      output_handle.write(pprint.pformat(sample_list))
    """
}


process rnaSPAdes {
    memory = 2.GB

    input:
      file reads from filteredPairedReads_ch2.collect()
      file datasets from rnaSPAdes_yaml_ch

    output:
      set val("rnaSPAdes"), \
          file("rnaSPAdes.gene_trans_map"), \
          file("rnaSPAdes/transcripts.fasta") \
          into rnaSPAdes_ch

    script:
      """
      rnaspades.py --dataset ${datasets} \
                   -t ${task.cpus} \
                   -m ${task.memory.toGiga()} \
                   -o rnaSPAdes \
                   --only-assembler -k 47 \
                   ${params.rnaSPAdes_params}

      # Make a fake gene to transcript file:
      cat "rnaSPAdes/transcripts.fasta" | grep ">" | tr -d ">" | cut -f 1 -d " " > tmp.txt
      paste tmp.txt tmp.txt > rnaSPAdes.gene_trans_map
      """
}

