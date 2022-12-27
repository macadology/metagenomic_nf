#!/usr/bin/env nextflow
// DSL 2 syntax
nextflow.enable.dsl=2

// help message
params.help = false
def helpMessage() {
    log.info"""
    ############################################################################
    ACRC kraken nextflow pipeline
    ----------------------------------------------------------------------------
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run kraken.nf

    Main input arguments:
        --querydir              Path to a folder containing all input fastq.
                                Required if readtype is (custom)
        --queryglob             Glob pattern of paired reads.
                                Required if readtype is (custom)
                                Optional for others.
                                E.g.: "*_{1,2}*{fastq,fastq.gz,fq,fq.gz}"
        --outputdir             Path to a folder containing all input fastq.
        --database              Database directory
        --overwrite             Overwrite existing directories. Otherwise,
                                skip. Default: false
        --profilers             fastp,decont,kraken,bracken


    Other input arguments:
    Look for editable params in nextflow.config, conf/*.conf and local.config where appropriate. Change the parameters either via the command line
    (e.g. --krakenDB), or by editing the config files.

    Decont:
    Removes human reads by mapping reads to GRCh38_latest_genomic.fna.gz and filtering out aligned reads.

    Default examples :
        # Show files detected
        nextflow run kraken.nf --querydir [Directory] --queryglob "*_{1,2}*{fastq,fastq.gz,fq,fq.gz}" --outputdir [Directory]

        # Run pipeline
        nextflow run kraken.nf -w [work directory] --querydir [Directory] --queryglob "*_{1,2}*{fastq,fastq.gz,fq,fq.gz}" --outputdir [Directory] --profilers fastp,decont,kraken2,bracken --database [path to kraken database]

    ############################################################################
    """.stripIndent()
}
if (params.help){
    helpMessage()
    exit 0
}

//=========== Load functions ============
// GroovyShell shell = new GroovyShell()
// def funcs = shell.parse(new File('./extra/functions.gvy'))
// (querydir, queryglob, outputdir) = funcs.parseio(params)

evaluate(new File("../extra/parseio.gvy"))

println ""
println "querydir : $querydir"
println "queryglob : $queryglob"
println "query : $querydir/**/$queryglob"
println "outputdir : $outputdir"
println "database: $params.database"
println ""

//============= Parse profilers ===========
Set profilers_expected = ['fastp', 'decont', 'kraken2', 'bracken']
Set profilers = []
params.profilers = ""
if(params.profilers.getClass() != Boolean){
  Set profilers_input = params.profilers.split(',')
  Set profiler_diff = profilers_input - profilers_expected
  profilers = profilers_input.intersect(profilers_expected)
  if( profiler_diff.size() != 0 ) {
  	log.warn "[Pipeline warning] Profiler $profiler_diff is not supported yet! Will only run $profilers.\n"
  }
}

println "Running softwares : $profilers"

//=========== Parameters ===========
include { FASTP } from '../modules/fastp'
include { DECONT } from '../modules/decontamination'
include { KRAKEN2; BRACKEN } from '../modules/kraken2'

workflow {
    //------ Get prefix ---------
    // Closure after Channel.fromFilePairs() is used to parse the prefix
    querydirname = file(querydir).name
    ch_input = Channel.fromFilePairs("$querydir/**/$queryglob", flat: false, size:params.size, maxDepth: params.maxdepth) { file ->
        for(int i = 0; i < 3 ; i++) {
            file = file.getParent()
            folder = file.getParent()
            // Look for prefix by going up the path structure
            if (querydirname == folder.name) {
                pref = file.name //prefix
                return pref
                break
            }
        }
    }

    //------- Display files ---------
    if(params.profilers.size() == 0 ) {
    	ch_input.view()
    }

    //------ Fastp + Decontamination -------
    if(profilers.contains('fastp')){
        FASTP(ch_input, outputdir)
        //FASTP.out.stdout.view { "FASTP STDOUT:\n$it" }
        ch_fastpreads = FASTP.out.reads
        //FASTP.out.reads.view()
    }else{
        ch_fastpreads = ch_input
    }

    if(profilers.contains('decont')){
        decontIndex = file(params.decontIndex)
        decontIndexDir = decontIndex.getParent()
        decontIndexName = decontIndex.getName()
        DECONT(ch_fastpreads, outputdir, decontIndexDir, decontIndexName)
        ch_reads = DECONT.out.reads
    }else{
        ch_reads = ch_fastpreads
    }

    //------ Kraken + Bracken ---------
    if(profilers.contains('kraken2')){
        KRAKEN2(ch_reads, outputdir, params.database)
    }

    if(profilers.contains('bracken')){
        if(profilers.contains('kraken2')){
            ch_kraken = KRAKEN2.out.output
        }else{
            dbname = file(params.database).name
            ch_kraken = Channel.fromFilePairs("$querydir/**/kraken2/$dbname/*.{kraken2.report,kraken2.tax}",flat: true, size: 2, checkIfExists: true) { file ->
                for(int i = 0; i < 3 ; i++) {
                    file = file.getParent()
                    folder = file.getParent()
                    // Look for prefix by going up the path structure
                    if (querydirname == folder.name) {
                        pref = file.name //prefix
                        return pref
                        break
                    }
                }
            }
            //ch_kraken.view()
        }
        BRACKEN(ch_kraken, outputdir, params.database)

    }

    // if(profilers.contains('test')){
    //     println "master: $outputdir"
    //     TEST(ch_reads, outputdir, params.testDatabase)
    //     TEST.out.stdout.view()
    // }
}
