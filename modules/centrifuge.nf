// run centrifuge
// Note, for the nt database, you need ~128gb ram to run, so it doesn't work properly locally for now.

process CENTRIFUGE {
    tag "${prefix}"
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }
    //container 'macadology/centrifuge'

    input:
    tuple val(prefix), path(reads1), path(reads2)
    path(centrifugeDBdir)

    output:
    publishDir "$params.outputdir/${prefix}/centrifuge", mode: 'copy'
    tuple path("${prefix}.centrifuge.summary.tsv"), path("${prefix}.centrifuge.kreport"), path("${prefix}.centrifuge.classification.out.tar.gz"), emit: output
    val("${prefix}"), emit: prefix
    stdout emit: stdout

    script:
    def outputdir = file("$params.procdir/${prefix}/centrifuge")
    if (outputdir.exists() && !params.overwrite) {
        println "$outputdir exists. Skipping $prefix ..."
        """
        exit 148
        """
    }else{
        """
        centrifuge -x $centrifugeDBdir/$params.centrifugeDBname -1 $reads1 -2 $reads2 -t -p ${task.cpus} -S ${prefix}.centrifuge.classification.out --report-file ${prefix}.centrifuge.summary.tsv

        centrifuge-kreport -x $centrifugeDBdir/$params.centrifugeDBname ${prefix}.centrifuge.classification.out > ${prefix}.centrifuge.kreport

        tar -zcvf ${prefix}.centrifuge.classification.out.tar.gz ${prefix}.centrifuge.classification.out

        rm ${prefix}.centrifuge.classification.out
        """
    }
}

// > ${prefix}.centrifuge.summary.tsv
// > ${prefix}.centrifuge.kreport
// > ${prefix}.centrifuge.classification.out.tar.gz
