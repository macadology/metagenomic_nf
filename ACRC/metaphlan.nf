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
    nextflow run metaphlan.nf

    Main input arguments:
        --querydir              Path to a folder containing all input fastq.
                                Required if readtype is (custom)
        --queryglob             Glob pattern of paired reads.
                                Required if readtype is (custom)
                                Optional for others.
                                E.g.: "*_{1,2}*{fastq,fastq.gz,fq,fq.gz}"
        --outputdir             Path to a folder containing all input fastq.
        --overwrite             Overwrite existing directories. Otherwise,
                                skip. Default: false
        --dryrun                List files detected. Default: false


    Other input arguments:
    Look for editable params in nextflow.config, conf/*.conf and local.config where appropriate. Change the parameters either via the command line
    (e.g. --krakenDB), or by editing the config files.

    Default examples :
        # Show files detected
        nextflow run metaphlan.nf --querydir [Directory] --queryglob "*_{1,2}*{fastq,fastq.gz,fq,fq.gz}" --outputdir [Directory] --dryrun true
        # Run program
        nextflow run metaphlan.nf --querydir [Directory] --queryglob "*_{1,2}*{fastq,fastq.gz,fq,fq.gz}" --outputdir [Directory]

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
println "outputdir : $outputdir"
println ""

//============= Parse profilers ===========
Set profilers_expected = ['fastp', 'decont', 'kraken2', 'bracken']
Set profilers = []
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
include { FASTP } from './modules/fastp'
include { DECONT } from './modules/decontamination'
include { KRAKEN2; BRACKEN } from './modules/kraken2'

workflow {
    //------ Get prefix ---------
    //// Closure after Channel.fromFilePairs() is used to parse the prefix
    //println "$querydir/**/$queryglob"
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
        //ch_reads = DECONT.out.reads
    }else{
        ch_reads = ch_fastpreads
    }

    //-------- Alignment ----------
    if(profilers.contains('BWA')){
        bwaIndex = file(params.bwaIndex)
        bwaIndexDir = bwaIndex.getParent()
        bwaIndexName = bwaIndex.getName()
        BWA(ch_reads, outputdir, bwaIndexDir, bwaIndexName)
        ch_sam = BWA.out.sam
        //MAPPED(BWA.out.output.collect(), outputdir, bwaIndexName)
        //ALIGN.out.stdout.view { "ALIGN STDOUT:\n$it" }
        //MAPPED.out.stdout.view { "MAPPED STDOUT:\n$it" }
    }else{
        ch_sam = ch_input
    }

    if(profilers.contains('bowtie')){
        btIndex = file(params.btIndex)
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

    if(profilers.contains('minimap')){
        // TBD: Add check database
        mmIndex = file(params.mmIndex)
        mmIndexDir = mmIndex.getParent()
        mmIndexName = mmIndex.getName()
        MINIMAP(ch_reads, outputdir, mmIndexDir, mmIndexName)
        ch_sam = MINIMAP.out.sam
        MINIMAP.out.stdout.view { "MINIMAP STDOUT:\n$it" }
    }else{
        ch_sam = ch_input
    }

    if(profilers.contains('sam2fastq')){
        SAM2FASTQ(ch_sam, outputdir)
        ch_reads = SAM2FASTQ.out.reads
    }

    //------ Kraken + Bracken ---------
    if(profilers.contains('kraken2')){
        if(params.krakenMMAP){
            KRAKEN2(ch_reads, outputdir, DB)
            // Note, this causes the output folder to be named 'database'. Fix it
        }else{
            KRAKEN2(ch_reads, outputdir, params.krakenDB)
        }
        //KRAKEN2.out.stdout.view { "KRAKEN2 STDOUT:\n$it" }
    }

    if(profilers.contains('bracken')){
        if(profilers.contains('kraken2')){
            ch_kraken = KRAKEN2.out.output
            //ch_kraken.view()
        }else{
            //ch_kraken = Channel.fromFilePairs("$params.procdir/**/kraken2/*.{report,tax}", flat: true, size: 2, checkIfExists: true, maxDepth: 2) { file -> file.getParent().getParent().name }
            dbname = file(params.krakenDB).name
            ch_kraken = Channel.fromFilePairs("$params.procdir/**/kraken2/$dbname/*.{kraken2.report,kraken2.tax}",flat: true, size: 2, checkIfExists: true) { file ->
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
        BRACKEN(ch_kraken, outputdir, params.krakenDB)
        //BRACKEN.out.stdout.view { "BRACKEN STDOUT:\n$it" }
    }

    //------ Centrifuge --------
    if(profilers.contains('centrifuge')){
        CENTRIFUGE(ch_reads, outputdir, params.centrifugeDBdir)
    }

    //------ Metaphlan3 + Humann3 -------
    if(profilers.contains('metaphlan')){
        METAPHLAN(ch_sam, outputdir, params.metaphlanDB_index, params.metaphlanDB_bt2Chocophlan)
        METAPHLAN.out.stdout.view()
    }

    if(profilers.contains('humann3')){
        HUMANN3(ch_reads, outputdir, params.humannDB_index, params.humannDB_Uniref, params.humannDB_Chocophlan, params.humannDB_bt2Chocophlan)
        HUMANN3.out.stdout.view()
    }

    //------ Test -------
    if(profilers.contains('test')){
        println "master: $outputdir"
        TEST(ch_reads, outputdir, params.testDatabase)
        TEST.out.stdout.view()
    }

    //------ srst2 -------
    if(profilers.contains('srst2')){
        SRST2(ch_reads, outputdir, params.srst2DB)
    }

    //------ Megares (RGI) ----------
    if(profilers.contains('rgi')){
        RGI(ch_reads, outputdir)
    }

    //------ diamond + megan
}