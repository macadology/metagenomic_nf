// run srst2

process SRST2 {
    tag "${prefix}"
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }
    container 'macadology/srst2'

    input:
    tuple val(prefix), path(reads1), path(reads2)
    path(srst2DB)

    output:
    publishDir "$params.procdir/${prefix}/srst2/${srst2DB.name}", mode: 'copy'
    path("${prefix}*"), emit: output
    val("${prefix}"), emit: prefix
    stdout emit: stdout

    script:
    def outputdir = file("$params.procdir/${prefix}/srst2/${srst2DB.name}")
    if (outputdir.exists() && !params.overwrite) {
        println "$outputdir exists. Skipping $prefix ..."
        """
        exit 148
        """
    }else{
        """
        srst2 --input_pe $reads1 $reads2 --output $prefix --log --gene_db $srst2DB --threads $params.srst2Threads
        """
    }
}

// > ${prefix}.log
// > ${prefix}_fastp.ARGannot_r3.pileup
// > ${prefix}_fastp.ARGannot_r3.sorted.bam
// > ${prefix}_fullgenes__ARGannot_r3__results.txt
// > ${prefix}_genes__ARGannot_r3__results.txt
