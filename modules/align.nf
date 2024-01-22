// run bwa mem
process BWA {
    // Align reads to target genome and count number of mapped reads
    tag "${prefix}"
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }
    //container 'macadology/bwa'

    input:
    tuple val(prefix), path(reads)
    val(procdir)
    path(bwaIndexDir)
    val(bwaIndexName)

    output:
    publishDir "$procdir/${prefix}/bwa", mode: 'copy'
    tuple val(prefix), path("${prefix}.${bwaIndexName}.sorted.bam"), emit: sam
    path "${prefix}.${bwaIndexName}.sorted.bam*", emit: output
    val("${prefix}"), emit: prefix
    stdout emit: stdout

    script:
    //String indexname = "${bwaIndex.baseName}_${bwaIndex.extension}"
    def outputdir = file("$procdir/${prefix}/bwa")
    if (outputdir.exists() && !params.overwrite) {
        println "$outputdir exists. Skipping $prefix ..."
        """
        exit 148
        """
    }else{
        """
        bwa mem -t ${task.cpus} ${bwaIndexDir}/${bwaIndexName} $reads | samtools sort - > ${prefix}.${bwaIndexName}.sorted.bam
        samtools index ${prefix}.${bwaIndexName}.sorted.bam > ${prefix}.${bwaIndexName}.sorted.bam.bai
        samtools view -c -F 268 ${prefix}.${bwaIndexName}.sorted.bam > ${prefix}.${bwaIndexName}.sorted.bam.mapped_counts
        #bwa mem -t ${task.cpus} ${bwaIndexDir}/${bwaIndexName} $reads | \\
        #samtools view $params.samtoolsOptions - > ${prefix}.${bwaIndexName}.sam
        #samtools view -c -F 268 ${prefix}.${bwaIndexName}.sam > ${prefix}.${bwaIndexName}.mapped_count
        """
    }
}

process BOWTIE {
    // Align reads to target genome and count number of mapped reads
    // 5th Dec, Set options......
    tag "${prefix}"
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }
    //container 'macadology/bwa'

    input:
    tuple val(prefix), path(reads)
    val(procdir)
    path(btIndexDir)
    val(btIndexName)

    output:
    publishDir "$procdir/${prefix}/bowtie", mode: 'copy'
    tuple val(prefix), path("${prefix}.${btIndexName}.sam"), emit: sam
    path "${prefix}.${btIndexName}.sam", emit: output
    val("${prefix}"), emit: prefix
    stdout emit: stdout

    script:
    //String indexname = "${btIndex.baseName}_${bwaIndex.extension}"
    reads1 = reads[0]
    reads2 = reads[1]
    def outputdir = file("$procdir/${prefix}/bowtie")
    if (outputdir.exists() && !params.overwrite) {
        println "$outputdir exists. Skipping $prefix ..."
        """
        exit 148
        """
    }else{
        """
        bowtie2 -p 12 --no-unal --no-head --no-sq --quiet --seed 1992 --very-sensitive -S ${prefix}.${btIndexName}.sam -x $btIndexDir/$btIndexName -1 $reads1 -2 $reads2 -S ${prefix}.${btIndexName}.sam
        """
    }
}

process MINIMAP {
    // Align reads to target genome and count number of mapped reads
    tag "${prefix}"
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }
    //container 'macadology/bwa'

    input:
    tuple val(prefix), path(reads)
    val(procdir)
    val(mmReadType)
    path(mmIndexDir)
    val(mmIndexName)

    output:
    publishDir "$procdir/${prefix}/minimap2", mode: 'copy'
    tuple val(prefix), path("${prefix}.${mmIndexName}.sam*"), emit: sam
    //path "${prefix}.${mmIndexName}.sorted.bam*", emit: output
    val("${prefix}"), emit: prefix
    stdout emit: stdout

    script:
    def outputdir = file("$procdir/${prefix}/minimap")
    if (outputdir.exists() && !params.overwrite) {
        println "$outputdir exists. Skipping $prefix ..."
        """
        exit 148
        """
    }else{
        //Note: Depending on whether minimap is used to align to a single fasta reference or multiple fasta reference, it'll affect the header output of the same file and affect subsequent bam conversion.
        """
        minimap2 -t ${task.cpus} -ax ${mmReadType} --split-prefix=tmp ${mmIndexDir}/${mmIndexName} $reads -o ${prefix}.${mmIndexName}.sam
        samtools view -c -F 268 ${prefix}.${mmIndexName}.sam > ${prefix}.${mmIndexName}.sam.mapped_counts
        #minimap2 -t ${task.cpus} -ax ${mmReadType} ${mmIndexDir}/${mmIndexName} $reads -o ${prefix}.${mmIndexName}.sam
        #samtools sort ${prefix}.${mmIndexName}.sam > ${prefix}.${mmIndexName}.sorted.bam
        #samtools index ${prefix}.${mmIndexName}.sorted.bam > ${prefix}.${mmIndexName}.sorted.bam.bai
        #samtools view -c -F 268 ${prefix}.${mmIndexName}.sorted.sam > ${prefix}.${mmIndexName}.sorted.bam.mapped_counts
        """
    }
}

process SAM2FASTQ {
    // Convert sam file to fastq files. Choose reads that align to ref.
    // 5th Dec 2022. Sam file requires header. Current code not working...
    // Only works for paired end reads for now....
    // Test if metaphlan works with sam file with header!
    tag "${prefix}"
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }

    input:
    tuple val(prefix), path(sam)
    val(procdir)

    output:
    publishDir "$procdir/${prefix}/sam2fastq/${refDB}", mode: 'copy'
    tuple val(prefix), path("${sam}_1.fq.gz"), path("${sam}_2.fq.gz"), emit: reads
    tuple path("${sam}_1.fq.gz"), path("${sam}_2.fq.gz"), emit: output
    val("${prefix}"), emit: prefix
    stdout emit: stdout

    script:
    def outputdir = file("$procdir/${prefix}/sam2fastq")
    if (outputdir.exists() && !params.overwrite) {
        println "$outputdir exists. Skipping $prefix ..."
        """
        exit 148
        """
    }else{
        """
        samtools fastq -F268 -1 ${sam}_1.fq -2 ${sam}_2.fq $sam
        gzip ${sam}_1.fq
        gzip ${sam}_2.fq
        """
    }

}

process MAPPED {
    // Create a csv file of counts of reads mapped to index
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }
    //container 'macadology/bwa'

    input:
    path(mapped_reads)
    val(procdir)
    val(bwaIndexName)

    output:
    publishDir "$procdir/reports", mode: 'copy'
    path "${bwaIndexName}.mappedreads.csv", emit: output
    stdout emit: stdout

    script:
    //String indexname = "${bwaIndex.baseName}_${bwaIndex.extension}"
    """
    > ${bwaIndexName}_mappedreads.csv
    for i in $mapped_reads; do echo \${i%%.*} \$(head \$i) >> ${bwaIndexName
        }.mappedreads.csv; done
    """
}
