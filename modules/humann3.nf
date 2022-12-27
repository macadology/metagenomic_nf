// run humann3
process METAPHLAN {
    tag "${prefix}"
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }
    //errorStrategy { task.exitStatus in 148 ? 'ignore' : 'ignore' }
    // http://huttenhower.sph.harvard.edu/humann_data/chocophlan/full_chocophlan.v201901_v31.tar.gz

    input:
    tuple val(prefix), path(reads)
    val(procdir)
    val(metaphlanDB_index)
    path(metaphlanDB_bt2Chocophlan)

    output:
    publishDir "$procdir/${prefix}/metaphlan", mode: 'move'
    tuple path("${prefix}_metaphlan_bugs_list.tsv"), emit: output
    val("${prefix}"), emit: prefix
    stdout emit: stdout

    script:
    // Note, existence check doesn't seem to work for aws right now...
    def outputdir = file("$procdir/${prefix}/metaphlan")
    //====================== Check sam vs fastq input ========================
    if (outputdir.exists() && !params.overwrite) {
        println "$outputdir exists. Skipping ${prefix} ..."
        """
        exit 148
        """
    }else{
        if (reads.size()==1) {
            // If input is a sam file or single ended read
            ext = (reads[0].getExtension()=="sam") ? "sam" : "fastq"
            """
            metaphlan ${reads} --bowtie2db $metaphlanDB_bt2Chocophlan --index $metaphlanDB_index --input_type $ext --nproc ${task.cpus} --min_alignment_len 100 -o ${prefix}_metaphlan_bugs_list.tsv
            """
        }else if (reads.size()==2) {
            // If onput is a paired end read
            """
            cat $reads > ${prefix}.fq.gz
            metaphlan ${prefix}.fq.gz --bowtie2db $metaphlanDB_bt2Chocophlan --index $metaphlanDB_index --input_type fastq --nproc ${task.cpus} --min_alignment_len 100 -o ${prefix}_metaphlan_bugs_list.tsv
            rm ${prefix}.fq.gz
            """
        }
    }
}

process HUMANN3 {
    tag "${prefix}"
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }
    //errorStrategy { task.exitStatus in 148 ? 'ignore' : 'ignore' }
    // http://huttenhower.sph.harvard.edu/humann_data/chocophlan/full_chocophlan.v201901_v31.tar.gz

    input:
    tuple val(prefix), path(reads)
    val(procdir)
    val(humannDB_index)
    path(humannDB_Uniref)
    path(humannDB_Chocophlan)
    path(humannDB_bt2Chocophlan)

    output:
    //publishDir "$procdir/${prefix}/humann3", mode: 'copy'
    publishDir "$procdir/${prefix}/humann3", mode: 'move' // Use only if no other processes use humann output files. Great for saving space on ACRC
    tuple val(prefix), path("${prefix}_*.tsv"), emit: output
    val("${prefix}"), emit: prefix
    stdout emit: stdout

    script:
    // Note, existence check doesn't seem to work for aws right now...
    def outputdir = file("$procdir/${prefix}/humann3")
    if (outputdir.exists() && !params.overwrite) {
        println "$outputdir exists. Skipping ${prefix} ..."
        """
        exit 148
        """
    }else{
        //$reads can be single-ended or paired-end. THe variable will expand automatically.
        """
        cat $reads > ${prefix}.fq.gz
        humann --input ${prefix}.fq.gz --output . --threads ${task.cpus} --protein-database $humannDB_Uniref --nucleotide-database $humannDB_Chocophlan --metaphlan-options '--bowtie2db $humannDB_bt2Chocophlan --index $humannDB_index --nproc ${task.cpus}' --bowtie-options '--threads ${task.cpus}' --diamond-options '--threads ${task.cpus}'
        rm ${prefix}.fq.gz
        mv ${prefix}*_humann_temp/${prefix}*_bowtie2_aligned.tsv .
        mv ${prefix}*_humann_temp/${prefix}*_diamond_aligned.tsv .
        mv ${prefix}*_humann_temp/${prefix}*.log .
        mv ${prefix}*_humann_temp/${prefix}*_metaphlan_bugs_list.tsv .
        rm -r ${prefix}*_humann_temp
        """
    }

}
