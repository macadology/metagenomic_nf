// run rgi

process RGI {
    tag "${prefix}"
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }
    //container 'macadology/rgi'

    input:
    tuple val(prefix), path(reads)
    val(procdir)

    output:
    publishDir "$procdir/${prefix}/rgi", mode: 'copy'
    path("${prefix}*"), emit: output
    val("${prefix}"), emit: prefix
    stdout emit: stdout

    script:
    def outputdir = file("$procdir/${prefix}/rgi")
    if (outputdir.exists() && !params.overwrite) {
        println "$outputdir exists. Skipping $prefix ..."
        """
        exit 148
        """
    }else{
        reads1 = reads[0]
        reads2 = reads[1]
        """
        which rgi
        rgi bwt -1 $reads1 -2 $reads2 -a bwa -n ${task.cpus} -o ${prefix}
        """
    }
}

// > ${prefix}.allele_mapping_data.json
// > ${prefix}.gene_mapping_data.txt
// > ${prefix}.sorted.temp.bam
// > ${prefix}.allele_mapping_data.txt
