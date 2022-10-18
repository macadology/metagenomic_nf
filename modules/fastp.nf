// decontamination using fastp

process FASTP {
    // validExitStatus 0,148
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }
    tag "${prefix}"

    input:
    tuple val(prefix), path(reads1), path(reads2)

    output:
    publishDir "$params.outputdir/${prefix}/fastp", mode: 'copy'
    tuple prefix, path("fastp_$reads1"), path("fastp_$reads2"), emit: reads
    tuple file("fastp_${prefix}.HTML"), file("fastp_${prefix}.json"), emit: output
    val("${prefix}"), emit: prefix
    stdout emit: stdout

    script:
    def outputdir = file("$params.outputdir/${prefix}/fastp")
    if (outputdir.exists() && !params.overwrite) {
        println "$outputdir exists. Skipping $prefix ..."
        """
        exit 148
        """
    }else{
        println "Executing fastp on $prefix ..."
        """
        fastp -q 25 -p 20 -i $reads1 -I $reads2 -o fastp_$reads1 -O fastp_$reads2 -j fastp_${prefix}.json -h fastp_${prefix}.HTML
        """
    }
    // > fastp_$reads1
    // > fastp_$reads2
    // > fastp_${prefix}.json
    // > fastp_${prefix}.HTML

}
