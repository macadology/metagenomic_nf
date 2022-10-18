// run bwa mem
process ALIGN {
    // Align reads to target genome and count number of mapped reads
    tag "${prefix}"
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }
    //container 'macadology/bwa'

    input:
    tuple val(prefix), path(reads1), path(reads2)
    path(bwaIndexDir)
    val(bwaIndexName)


    output:
    publishDir "$params.outputdir/${prefix}/bwa", mode: 'copy'
    path "${prefix}.${bwaIndexName}.out", emit: output
    val("${prefix}"), emit: prefix
    stdout emit: stdout

    script:
    //String indexname = "${bwaIndex.baseName}_${bwaIndex.extension}"
    def outputdir = file("$params.outputdir/${prefix}/bwa")
    if (outputdir.exists() && !params.overwrite) {
        println "$outputdir exists. Skipping $prefix ..."
        """
        exit 148
        """
    }else{
        """
        bwa mem -t ${task.cpus} ${bwaIndexDir}/${bwaIndexName} $reads1 $reads2 | \\
        samtools view $params.samtoolsOptions - > "${prefix}.${bwaIndexName}.out"
        #> ${prefix}.${bwaIndexName}.out
        """
    }
}

process MAPPED {
    // Create a csv file of counts of reads mapped to index
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }
    //container 'macadology/bwa'

    input:
    path(mapped_reads)
    val(bwaIndexName)

    output:
    publishDir "$params.outputdir/reports", mode: 'copy'
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
