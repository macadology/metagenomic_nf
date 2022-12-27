 // decontamination
process DECONT {
    // Remove reads that are aligned to target sequence
    // TBD: Only works for paired reads for now. Add long read decont
    tag "${prefix}"
    errorStrategy { task.exitStatus in 148 ? 'ignore' : params.errorstrategy }

    input:
    tuple val(prefix), path(reads)
    val(procdir)
    path(decontIndexDir)
    val(decontIndexName)

    output:
    publishDir "$procdir/${prefix}/fastp", mode: 'copy'
    tuple val(prefix), path("decont_*"), emit: reads
    path("aln-se_${prefix}.sam")
    val("${prefix}"), emit: prefix
    stdout emit: stdout

    script:
    if (reads.size()==2) {
        // For paired end reads
        reads1 = reads[0]
        reads2 = reads[1]
        decontread1 = new File("$params.procdir/${prefix}/fastp/decont_${decontIndexName}_$reads1")
        decontread2 = new File("$params.procdir/${prefix}/fastp/decont_${decontIndexName}_$reads2")
        println reads1
        println reads2

        if (decontread1.exists() && decontread2.exists() && !params.overwrite) {
            println "$decontread1 and $decontread2 exists. Skipping $prefix ..."
            """
            exit 148
            """
        }else{
            """
            bwa mem -t ${task.cpus} $decontIndexDir/${decontIndexName} $reads1 $reads2 > aln-se_${prefix}.sam
            samtools fastq -f12 -F256  -1  decont_${decontIndexName}_$reads1 -2 decont_${decontIndexName}_$reads2 aln-se_${prefix}.sam
            """
        }
    } else if (reads.size()==1) {
        // For single read (long reads)
        decontread = new File("$params.procdir/${prefix}/fastp/decont_${decontIndexName}_$reads")
        if (decontread.exists() && !params.overwrite) {
            println "$decontread exist. Skipping $prefix ..."
            """
            exit 148
            """
        }else{
            """
            bwa mem -t ${task.cpus} $decontIndexDir/${decontIndexName} $reads > aln-se_${prefix}.sam
            samtools fastq -f12 -F256 decont_${decontIndexName}_$reads aln-se_${prefix}.sam
            """
        }
    }

}
