process TEST {
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }
    container 'ubuntu:18.04'

    input:
    tuple val(prefix), path(reads1), path(reads2)
    path(testDatabase)

    output:
    publishDir "$params.procdir/${prefix}/test", mode: 'copy'
    path("${prefix}.txt"), emit: output
    stdout emit: stdout
    val("${prefix}"), emit: prefix

    script:
    def outputdir = file("$params.procdir/${prefix}/test")
    if (outputdir.exists()) {
        println "$outputdir exists. Overwriting $prefix ..."
    }
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
