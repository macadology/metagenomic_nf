process GEN_SR {
    tag "${prefix}"
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }
    //errorStrategy { task.exitStatus in 148 ? 'ignore' : 'ignore' }
    // http://huttenhower.sph.harvard.edu/humann_data/chocophlan/full_chocophlan.v201901_v31.tar.gz

    input:
    tuple val(prefix), path(reads)
    val(procdir)

    output:
    publishDir "$procdir/${prefix}/metaphlan", mode: 'copy'
    tuple val(prefix), path("${prefix}*.fastq"), emit: output
    path("${prefix}_LR2SR.fasta"), emit: fasta
    val("${prefix}"), emit: prefix
    stdout emit: stdout

    script:
    // Note, existence check doesn't seem to work for aws right now...
    def outputdir = file("$procdir/${prefix}/wgsim")
    //====================== Check sam vs fastq input ========================
    if (outputdir.exists() && !params.overwrite) {
        println "$outputdir exists. Skipping ${prefix} ..."
        """
        exit 148
        """
    }else{
        println "Generating Short Reads ${prefix} ..."
        """
        seqtk seq -a ${reads} > ${prefix}_LR2SR.fasta
        counts=\$(grep -v ">" *.fasta | wc -c)
        Nreads=\$((\$counts/100*10))
        echo "Nreads=\$Nreads"
        wgsim -e 0 -N \$Nreads -d 1 -s 1 -1 100 -2 100 -r 0 -R 0.0001 -X 0.0001 -S 123 ${prefix}_LR2SR.fasta ${prefix}_LR2SR_R1.fastq ${prefix}_LR2SR_R2.fastq
        """

    }
}
