process TEST {
    tag "${prefix}"
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }

    input:
    tuple val(prefix), file(reads)
    val(procdir)
    val(testDatabase)

    output:
    //publishDir "$procdir/${prefix}/test", mode: 'copy'
    //path("${prefix}.txt"), emit: output
    stdout emit: stdout
    val("${prefix}"), emit: prefix

    script:
    outputdir = file("$procdir/${prefix}/test")
    reads1 = reads[0]
    reads2 = reads[1]
    """
    echo $prefix $reads
    """
    println "variable type: ${reads.getClass()}"
    println "array size: ${reads.size()}"
    println "read1 suffix : ${reads[0].getExtension()}"
    // Note
    //file('random') returns ~/random
    //file('/random') returns /random
    if (outputdir.exists() && !params.overwrite) {
        println "$outputdir exists. Overwriting $prefix ..."
        """
        exit 148
        """
    }else{
        """
        echo "${prefix} sleeping for 10 sec."
        #sleep 10
        #cat $testDatabase
        #touch ${prefix}.txt
        #head -n 1 $reads >> ${prefix}.txt
        #head -n 1 $reads >> ${prefix}.txt
        #touch ${prefix}.txt
        """
    }
}
