// run bwa mem
process BWA {
    // Align reads to target genome and count number of mapped reads
    tag "${prefix}"
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }
    //container 'macadology/bwa'

    input:
    tuple val(prefix), path(reads1), path(reads2)
    val(procdir)
    path(bwaIndexDir)
    val(bwaIndexName)

    output:
    publishDir "$procdir/${prefix}/bwa", mode: 'copy'
    tuple val(prefix), path("${prefix}.${bwaIndexName}.sam"), emit: sam
    path "${prefix}.${bwaIndexName}.sam", emit: output
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
        bwa mem -t ${task.cpus} ${bwaIndexDir}/${bwaIndexName} $reads1 $reads2 | \\
        samtools view $params.samtoolsOptions - > "${outputdir}/${prefix}.${bwaIndexName}.sam"
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
    tuple val(prefix), path(reads1), path(reads2)
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
    tuple val(prefix), path(reads1), path(reads2)
    val(procdir)
    path(mmIndexDir)
    val(mmIndexName)

    output:
    publishDir "$procdir/${prefix}/minimap2", mode: 'copy'
    tuple val(prefix), path("${prefix}.${mmIndexName}.sam"), emit: sam
    path "${prefix}.${mmIndexName}.sam", emit: output
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
        """
        minimap2 -t ${task.cpus} -ax sr ${mmIndexDir}/${mmIndexName} ${reads1} ${reads2} > ${prefix}.${mmIndexName}.sam
        """
    }
}

process MINIMAP_LR {
    // Align reads to target genome and count number of mapped reads
    tag "${prefix}"
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }
    //container 'macadology/bwa'

    input:
    tuple val(prefix), path(read)
    val(procdir)
    path(mmIndexDir)
    val(mmIndexName)

    output:
    publishDir "$procdir/${prefix}/minimap2", mode: 'copy'
    tuple val(prefix), path("${prefix}.${mmIndexName}.sam"), emit: sam
    path "${prefix}.${mmIndexName}.sam", emit: output
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
        """
        minimap2 -t ${task.cpus} -ax ${params.mmreadtype} ${mmIndexDir}/${mmIndexName} ${read} > ${prefix}.${mmIndexName}.sam
        """
    }
}

process SAM2FASTQ {
    // Convert sam file to fastq files. Choose reads that align to ref.
    // 5th Dec 2022. Sam file requires header. Current code not working...
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

// PREFIX=TR
//for i in $(find . -type f -name "*.$PREFIX.fasta.txt"); do echo $(basename $(dirname $(dirname $i))) $(head $i) >> ~/GIS/food_fermentation/reports/$PREFIX.fasta.mappedreads.csv; done
//for i in $(find . -type f -name "*.AH.fasta.txt"); do echo $(head $i); done
