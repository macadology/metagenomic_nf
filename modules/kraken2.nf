// run kraken and bracken
// https://www.researchgate.net/profile/Juan-Jovel/publication/338686214_How_to_run_Kraken_on_metagenomics_FASTQ_files/links/5e24b527a6fdcc1015781691/How-to-run-Kraken-on-metagenomics-FASTQ-files.pdf
// https://hpc.nih.gov/apps/kraken.html
// https://github.com/DerrickWood/kraken2/issues/120

mmap = ""
if(params.krakenMMAP){
    mmap = "--memory-mapping "
}

process KRAKEN2 {
    tag "${prefix}"
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }
    container 'macadology/kraken'

    input:
    tuple val(prefix), path(reads)
    val(procdir)
    path(krakenDB)

    output:
    publishDir "$procdir/${prefix}/kraken2/${krakenDB.name}", mode: 'copy'
    tuple val("${prefix}"), path("${prefix}.kraken2.report"), path("${prefix}.kraken2.tax"), emit: output
    val("${prefix}"), emit: prefix
    stdout emit: stdout

    script:
    // For paired end reads, $paired will add the --paired option and $reads will expand automatically.
    def paired = (reads.size()==2) ? "--paired " : ""
    def outputdir = file("$procdir/${prefix}/kraken2/${krakenDB.name}")
    if (outputdir.exists() && !params.overwrite) {
        println "$outputdir exists. Skipping $prefix ..."
        """
        exit 148
        """
    }else{
        """
        kraken2 $mmap--db $krakenDB $paired--threads $params.krakenThreads --output ${prefix}.kraken2.out --report ${prefix}.kraken2.tax $reads --use-mpa-style > /dev/null

        kraken2 $mmap--db $krakenDB $paired--threads $params.krakenThreads --report ${prefix}.kraken2.report $reads > /dev/null

        if [ ! ${params.krakenKeepOutput} == "true" ]; then
            rm ${prefix}.kraken2.out
        fi
        """
    }
}

process BRACKEN {
    tag "${prefix}"
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }
    container 'macadology/kraken2'

    input:
    tuple val(prefix), path(krakenReport), path(krakenTax)
    val(procdir)
    path(brackenDB)

    output:
    publishDir "$procdir/${prefix}/kraken2/${brackenDB.name}", mode: 'copy'
    tuple path("${prefix}.bracken.P"), path("${prefix}.bracken.F"), path("${prefix}.bracken.G"), path("${prefix}.bracken.S"), emit: output
    stdout emit: stdout
    val("${prefix}"), emit: prefix

    script:
    def outputdir = file("$procdir/${prefix}/kraken2/${brackenDB.name}/${prefix}.bracken.P")
    if (outputdir.exists() && !params.overwrite) {
        println "$outputdir exists. Skipping $prefix ..."
        """
        exit 148
        """
    }else{
        """
        which bracken
        for taxa in P F G S
        do
            bracken -d $brackenDB -r $params.krakenReadlength -i $krakenReport -o ${prefix}.bracken.\$taxa -l \$taxa
        done
        """
    }
}
//docker run $(for i in $(ls $PWD); do echo " -v $(readlink -f $i):/home/ubuntu/$(basename $i)"; done) -it macadology/kraken /bin/bash
