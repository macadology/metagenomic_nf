// run kraken and bracken
// todo, store database in ram, run all kraken classification first, then run subsequent analysis
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
    tuple prefix, path(reads1), path(reads2)
    path(krakenDB)

    output:
    publishDir "$params.procdir/${prefix}/kraken2/${krakenDB.name}", mode: 'copy'
    tuple val("${prefix}"), path("${prefix}.kraken2.report"), path("${prefix}.kraken2.tax"), emit: output
    val("${prefix}"), emit: prefix
    stdout emit: stdout

    script:
    def outputdir = new File("$params.procdir/${prefix}/kraken2/${krakenDB.name}")
    if (outputdir.exists() && !params.overwrite) {
        println "$outputdir exists. Skipping $prefix ..."
        """
        exit 148
        """
    }else{
        """
        kraken2 $mmap--db $krakenDB --paired --threads $params.krakenThreads --output ${prefix}.kraken2.out --report ${prefix}.kraken2.tax $reads1 $reads2 --use-mpa-style > /dev/null

        kraken2 $mmap--db $krakenDB --paired --threads $params.krakenThreads --report ${prefix}.kraken2.report $reads1 $reads2 > /dev/null

        if [ ! ${params.krakenKeepOutput} == "true" ]; then
            rm ${prefix}.kraken2.out
        fi
        #> ${prefix}.kraken2.report
        #> ${prefix}.kraken2.tax
        """
    }
}

process BRACKEN {
    tag "${prefix}"
    container 'macadology/kraken2'

    input:
    tuple val(prefix), path(krakenReport), path(krakenTax)
    path(brackenDB)

    output:
    publishDir "$params.procdir/${prefix}/kraken2/${brackenDB.name}", mode: 'copy'
    tuple path("${prefix}.bracken.P"), path("${prefix}.bracken.F"), path("${prefix}.bracken.G"), path("${prefix}.bracken.S"), emit: output
    stdout emit: stdout
    val("${prefix}"), emit: prefix

    script:
    """
    which bracken
    for taxa in P F G S
    do
        bracken -d $brackenDB -r $params.krakenReadlength -i $krakenReport -o ${prefix}.bracken.\$taxa -l \$taxa
    done
    """
}
//docker run $(for i in $(ls $PWD); do echo " -v $(readlink -f $i):/home/ubuntu/$(basename $i)"; done) -it macadology/kraken /bin/bash
