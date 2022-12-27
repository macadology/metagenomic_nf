// decontamination using fastp

process FASTP {
    // validExitStatus 0,148
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }
    tag "${prefix}"

    input:
    tuple val(prefix), path(reads)
    val(procdir)

    output:
    publishDir "$procdir/${prefix}/fastp", mode: 'copy'
    tuple val(prefix), path("fastp_$reads1"), path("fastp_$reads2"), emit: reads
    tuple file("fastp_${prefix}.HTML"), file("fastp_${prefix}.json"), emit: output
    val("${prefix}"), emit: prefix
    stdout emit: stdout

    script:
    def outputdir = file("$procdir/${prefix}/fastp")
    if (outputdir.exists() && !params.overwrite) {
        println "$outputdir exists. Skipping $prefix ..."
        """
        exit 148
        """
    }else{
        println "Executing fastp on $prefix ..."
        if (reads.size()==2) {
            reads1 = reads[0]
            reads2 = reads[1]
            """
            fastp -q 25 -p 20 -i $reads1 -I $reads2 -o fastp_$reads1 -O fastp_$reads2 -j fastp_${prefix}.json -h fastp_${prefix}.HTML
            """
        } else (reads.size()==1) {
            """
            fastp -q 25 -p 20 -i $reads -o fastp_$reads1 -O fastp_$reads2 -j fastp_${prefix}.json -h fastp_${prefix}.HTML
            """
        }
    }
}
