// run rgi

process RGI {
    tag "${prefix}"
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }
    container 'macadology/rgi'

    input:
    tuple val(prefix), path(reads1), path(reads2)

    output:
    publishDir "$params.procdir/${prefix}/rgi", mode: 'copy'
    path("${prefix}*"), emit: output
    val("${prefix}"), emit: prefix
    stdout emit: stdout

    script:
    def outputdir = file("$params.procdir/${prefix}/rgi")
    if (outputdir.exists() && !params.overwrite) {
        println "$outputdir exists. Skipping $prefix ..."
        """
        exit 148
        """
    }else{
        """
        which rgi
        rgi bwt -1 $reads1 -2 $reads2 -a bwa -n $params.rgiThreads -o ${prefix}
        """
    }
}

// > ${prefix}.allele_mapping_data.json
// > ${prefix}.gene_mapping_data.txt
// > ${prefix}.sorted.temp.bam
// > ${prefix}.allele_mapping_data.txt
