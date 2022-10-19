process TEST {
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }

    input:
    tuple val(prefix), path(reads1), path(reads2)
    val(procdir)
    //path(procdir) doesn't work because procdir is not a file
    path(testDatabase)

    output:
    publishDir "$procdir/${prefix}/test", mode: 'copy'
    path("${prefix}.txt"), emit: output
    stdout emit: stdout
    val("${prefix}"), emit: prefix

    script:
    outputdir = file("$procdir/${prefix}/test")
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
        cat $testDatabase
        touch ${prefix}.txt
        head -n 1 $reads1 >> ${prefix}.txt
        head -n 1 $reads2 >> ${prefix}.txt
        touch ${prefix}.txt
        """
    }
}
