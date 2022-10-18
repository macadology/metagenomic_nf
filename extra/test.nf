#!/usr/bin/env nextflow
// DSL 2 syntax
nextflow.enable.dsl=2

// help message
params.help = false
def helpMessage() {
    log.info"""
    Testing Aws batch
    """.stripIndent()
}

//---------- Define process -----------
process PRINTHEAD {
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }
    container 'ubuntu:18.04'
    queue 'jon-bioinfo-queue'

    input:
    tuple val(prefix), path(reads1), path(reads2)

    output:
    publishDir "$params.procdir/${prefix}/test", mode: 'copy'
    path("${prefix}.txt"), emit: output
    stdout emit: stdout
    val("${prefix}"), emit: prefix

    script:
    def outputdir = new File("$params.procdir/${prefix}/test")
    if (outputdir.exists()) {
        println "$outputdir exists. Skipping $prefix ..."
        """
        exit 148
        """
    }else{
        """
        echo ${prefix}
        > ${prefix}.txt
        head -n 1 $reads1 >> ${prefix}.txt
        head -n 1 $reads2 >> ${prefix}.txt
        head -n 1 $reads1
        head -n 1 $reads2
        """
    }
}

//---------- Main workflow ------------
workflow {
    querydir = params.procdir
    queryglob = "fastp_*_{1,2}*{fastq,fastq.gz,fq,fq.gz}"

    //------ Get Prefix ---------
    //println "$querydir/**/$queryglob"
    ch_input = Channel.fromFilePairs("$querydir/**/$queryglob", flat: true, maxDepth: 3) { file ->
        for(int i = 0; i < 2 ; i++) {
            file = file.getParent()
            pref = file.name //prefix
            if (pref ==~ /[A-Z]{3}\d{3}/){ //The pattern means 3 letters 3 nums.
                return pref
                break
            }
        }
    }
    //ch_input.view()
    PRINTHEAD(ch_input)
    //PRINTHEAD.out.stdout.view { "PRINTHEAD STDOUT:\n$it" }
}
