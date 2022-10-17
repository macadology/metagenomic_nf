 // decontamination

process DECONT {
    // Remove reads that are aligned to target sequence
    tag "${prefix}"
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }

    input:
    tuple val(prefix), path(reads1), path(reads2)
    path(decontIndexDir)
    val(decontIndexName)

    output:
    publishDir "$params.procdir/${prefix}/fastp", mode: 'copy'
    tuple val(prefix), path("decont_${decontIndexName}_$reads1"), path("decont_${decontIndexName}_$reads2"), emit: reads
    path("aln-se_${prefix}.sam")
    val("${prefix}"), emit: prefix
    stdout emit: stdout

    script:
    //def matcher = "$params.decontIndex" =~ /([^\/]+)\.(.*)$/
    //println "${matcher[0][1]}"
    //String indexname = "${decontIndex.baseName}_${decontIndex.extension}"
    def decontread1 = new File("$params.procdir/${prefix}/fastp/decont_${decontIndexName}_$reads1")
    def decontread2 = new File("$params.procdir/${prefix}/fastp/decont_${decontIndexName}_$reads2")
    if (decontread1.exists() && decontread2.exists() && !params.overwrite) {
        println "$decontread1 and $decontread2 exists. Skipping $prefix ..."
        """
        exit 148
        """
    }else if (params.test) {
        println "(Test) Running decont on ${prefix}"
        """
        exit 148
        """
    }else{
        """
        #which bwa
        bwa mem -t $params.bwaThreads $decontIndexDir/${decontIndexName} $reads1 $reads2 > aln-se_${prefix}.sam
        samtools fastq -f12 -F256  -1  decont_${decontIndexName}_$reads1 -2 decont_${decontIndexName}_$reads2 aln-se_${prefix}.sam
        #> aln-se_${prefix}.sam
        #> decont_${decontIndexName}_$reads1
        #> decont_${decontIndexName}_$reads2
        """
    }

}
