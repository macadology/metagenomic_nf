// run humann3

process HUMANN3 {
    tag "${prefix}"
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }
    //errorStrategy { task.exitStatus in 148 ? 'ignore' : 'ignore' }
    //container 'biobakery/humann:3.1.1'
    // http://huttenhower.sph.harvard.edu/humann_data/chocophlan/full_chocophlan.v201901_v31.tar.gz
    // container 'macadology/humann3 '

    input:
    tuple val(prefix), path(reads1), path(reads2)
    val(procdir)
    path(humannDB_Uniref)
    path(humannDB_Chocophlan)
    path(humannDB_bt2Chocophlan)

    output:
    publishDir "$procdir/${prefix}/humann3", mode: 'copy'
    tuple path("${prefix}_*.tsv"), emit: output
    //, path("${prefix}*_humann_temp/${prefix}*")
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
    """
    #exit 148 #For testing purposes, please remove after test.
    which humann
    which metaphlan
    echo ${prefix}
    cat $reads1 $reads2 > ${prefix}.fq.gz
    humann --input ${prefix}.fq.gz --output . --threads ${task.cpus} --protein-database $humannDB_Uniref --nucleotide-database $humannDB_Chocophlan --metaphlan-options '--bowtie2db $humannDB_bt2Chocophlan --index $params.humannDB_index --nproc ${task.cpus}' --bowtie-options '--threads ${task.cpus}' --diamond-options '--threads ${task.cpus}'
    rm ${prefix}.fq.gz
    mv ${prefix}*_humann_temp/${prefix}_bowtie2_aligned.tsv .
    mv ${prefix}*_humann_temp/${prefix}_diamond_aligned.tsv .
    mv ${prefix}*_humann_temp/${prefix}.log .
    mv ${prefix}*_humann_temp/${prefix}_metaphlan_bugs_list.tsv .
    rm -r ${prefix}*_humann_temp
    """
    }

}

// > ${prefix}_genefamilies.tsv
// > ${prefix}_pathabundance.tsv
// > ${prefix}_pathcoverage.tsv
// mkdir ${prefix}_humann_temp
// > ${prefix}_humann_temp/${prefix}_bowtie2_aligned.sam
// > ${prefix}_humann_temp/${prefix}_bowtie2_aligned.tsv
