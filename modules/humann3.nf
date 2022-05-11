// run humann3

process HUMANN3 {
    tag "${prefix}"
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }
    container 'macadology/humann3'

    input:
    tuple prefix, path(reads1), path(reads2)
    path(humannDB_Uniref)
    path(humannDB_Chocophlan)
    path(humannDB_bt2Chocophlan)

    output:
    publishDir "$params.procdir/${prefix}/humann3", mode: 'copy'
    tuple path("${prefix}_*.tsv"), path("${prefix}*_humann_temp/${prefix}*"), emit: output
    val("${prefix}"), emit: prefix
    stdout emit: stdout

    script:
    def outputdir = new File("$params.procdir/${prefix}/humann3")
    if (outputdir.exists() && !params.overwrite) {
        println "$outputdir exists. Skipping $params ..."
        """
        exit 148
        """
    }else{
    """
    which humann
    which metaphlan
    cat $reads1 $reads2 > ${prefix}.fq.gz
    humann --input ${prefix}.fq.gz --output . --threads $params.humannThreads --protein-database $humannDB_Uniref --nucleotide-database $humannDB_Chocophlan --metaphlan-options '--bowtie2db $humannDB_bt2Chocophlan --nproc $params.humannThreads'
    rm ${prefix}.fq.gz
    """
    }

}

// > ${prefix}_genefamilies.tsv
// > ${prefix}_pathabundance.tsv
// > ${prefix}_pathcoverage.tsv
// mkdir ${prefix}_humann_temp
// > ${prefix}_humann_temp/${prefix}_bowtie2_aligned.sam
// > ${prefix}_humann_temp/${prefix}_bowtie2_aligned.tsv
