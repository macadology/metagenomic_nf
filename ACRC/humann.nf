#!/usr/bin/env nextflow
// DSL 2 syntax
nextflow.enable.dsl=2

// help message
params.help = false
def helpMessage() {
    log.info"""
    ############################################################################
    ACRC metaphlan nextflow pipeline
    ----------------------------------------------------------------------------
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run humann.nf

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
        --profilers             fastp,decont,metaphlan,humann


    Other input arguments:
    Look for editable params in nextflow.config, conf/*.conf and local.config where appropriate. Change the parameters either via the command line
    (e.g. --humannDB), or by editing the config files.

    Default examples :
        # Show files detected
        nextflow run humann.nf --querydir [Directory] --queryglob "*_{1,2}*{fastq,fastq.gz,fq,fq.gz}" --outputdir [Directory]

        # Run program
        nextflow run humann.nf --querydir [Directory] --queryglob "*_{1,2}*{fastq,fastq.gz,fq,fq.gz}" --outputdir [Directory] --profilers fastp,decont,humann --database [path to humann database]

    ############################################################################
    """.stripIndent()
}
if (params.help){
    helpMessage()
    exit 0
}

//=========== Input/Output ===========
// Assumes a particular folder structure. See README.md for more information
// Define (querydir) Directory containing reads and (queryglob) Read file format 

//------ queryglob, querydir, outputdir
if(params.outputdir == "" || params.querydir == ""){
    log.error "Error: You are running in general mode. Please specify a query (--querydir) and output (--outputdir) directory"
    exit 0
}

if(params.queryglob){
    queryglob = params.queryglob
}else{
    queryglob = "*_{1,2}*{fastq,fastq.gz,fq,fq.gz}"
}
querydir = params.querydir
outputdir = params.outputdir

println ""
println "querydir : $querydir"
println "queryglob : $queryglob"
println "query : $querydir/**/$queryglob"
println "outputdir : $outputdir"
println ""

//============= Parse profilers ===========
Set profilers_expected = ['fastp', 'bowtie', 'decont', 'metaphlan', 'humann']
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
include { BOWTIE } from '../modules/align'
include { METAPHLAN; HUMANN3 } from '../modules/humann3'

workflow {
    //------ Get prefix ---------
    // Closure after Channel.fromFilePairs() is used to parse the prefix
    querydirname = file(querydir).name
    ch_input = Channel.fromFilePairs("$querydir/**/$queryglob", flat: true, size:params.size, maxDepth: params.maxdepth) { file ->
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

    //------ Metaphlan3 + Humann3 -------
    if(profilers.contains('bowtie')){
        if((profilers.contains('metaphlan') || profilers.contains('humann')) && !params.btIndex){
            btIndex = file("${params.humannDB_bt2Chocophlan}/${params.humannDB_index}")
        } else {
            btIndex = file(params.btIndex)
        }
        btIndexDir = btIndex.getParent()
        btIndexName = btIndex.getName()
        BOWTIE(ch_reads, outputdir, btIndexDir, btIndexName)
        ch_sam = BOWTIE.out.sam
        //MAPPED(ALIGN.out.output.collect(), outputdir, bwaIndexName)
        //ALIGN.out.stdout.view { "ALIGN STDOUT:\n$it" }
        //MAPPED.out.stdout.view { "MAPPED STDOUT:\n$it" }
    }else{
        ch_sam = ch_input
    }

    if(profilers.contains('metaphlan')){
        METAPHLAN(ch_sam, outputdir, params.metaphlanDB_index, params.metaphlanDB_bt2Chocophlan)
        METAPHLAN.out.stdout.view()
    }

    if(profilers.contains('humann')){
        HUMANN3(ch_sam, outputdir, params.humannDB_index, params.humannDB_Uniref, params.humannDB_Chocophlan, params.humannDB_bt2Chocophlan)
        HUMANN3.out.stdout.view()
    }

    if(profilers.contains('test')){
        println "master: $outputdir"
        TEST(ch_reads, outputdir, params.testDatabase)
        TEST.out.stdout.view()
    }
}
