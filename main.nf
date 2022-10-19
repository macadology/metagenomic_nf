#!/usr/bin/env nextflow
// DSL 2 syntax
nextflow.enable.dsl=2

// help message
params.help = false
def helpMessage() {
    log.info"""
    ############################################################################
    Jon's metagenomic nextflow pipeline
    ----------------------------------------------------------------------------
    Refer to README.md for more details.
    Usage:
    The typical command for running the pipeline is as follows:
        nextflow run main.nf --profilers [software]

    Main input arguments:
        --general               Determines if pipeline runs with general
                                parameters (true) or project specific parameters in local.config (false). Default: false
        --querydir              Path to a folder containing all input fastq.
                                Required if readtype is (custom)
        --queryglob             Glob pattern of paired reads.
                                Required if readtype is (custom)
                                Optional for others.
                                E.g.: "*_{1,2}*{fastq,fastq.gz,fq,fq.gz}"
        --profilers             Metagenomcis profilers to run. If not defined,
                                show matched reads. Options -> 'fastp', 'decont', 'kraken2', 'bracken', 'centrifuge', 'humann3', 'srst2', 'rgi'
        --readtype              (raw) raw reads
                                (fastp) fastp reads *DEFAULT*
                                (decont) human decontaminated fastp reads
        --overwrite             Overwrite existing directories. Otherwise,
                                skip. Default: false
        --outputdir             Path to a folder containing all input fastq.

    Other input arguments:
    Look for editable params in nextflow.config, conf/*.conf and local.config where appropriate. Change the parameters either via the command line
    (e.g. --krakenDB), or by editing the config files.

    Default examples :
        nextflow run main.nf --general true --querydir [Directory] --queryglob "*_{1,2}*{fastq,fastq.gz,fq,fq.gz}" --outputdir [Directory] -profilers fastp,kraken2,bracken #General mode

    Jonai examples :
        nextflow run main.nf -profile jonai #Show the reads that'll be used as inputs
        nextflow run main.nf -profile jonai --readtype raw --profilers fastp,decont
        nextflow run main.nf -profile jonai -w [WORKDIR] --krakenMMAP --profilers kraken2,bracken #Preload kraken database. See README.md
        nextflow run main.nf -profile jonai --profilers humann3

    AWSBatch examples :
        nextflow run main.nf -profile batch --bucket-dir s3://jon-nextflow-work --profilers humann3

    ACRC examples :
        nextflow run main.nf -profile acrc --general true --querydir [Directory] --queryglob "*_{1,2}*{fastq,fastq.gz,fq,fq.gz}" --outputdir [Directory] -w [WORKDIR] # Show the reads that'll be used as inputs
        nextflow run main.nf -profile acrc --general true --querydir [Directory] --queryglob "*_{1,2}*{fastq,fastq.gz,fq,fq.gz}" --outputdir [Directory] -w [WORKDIR] --profilers humann3

    ############################################################################
    """.stripIndent()
}
if (params.help){
    helpMessage()
    exit 0
}

//=========== Input/Output ===========
// Assumes a particular folder structure. See README.md for more information
// Possible inputs --> (raw) raw reads, (fastp) fastp reads, (decont) human decontaminated fastp reads (custom) custom reads
// Define (querydir) Directory containing reads and (queryglob) Read file format 

if(params.general){
    //Checks if local.config is a placeholder or specific to project.
    if(params.outputdir == "" || params.querydir == ""){
        log.error "Error: You are running in general mode. Please specify a query (--querydir) and output (--outputdir) directory"
        exit 0
    }
    //------ queryglob, querydir, outputdir
    if(params.queryglob){
        queryglob = params.queryglob
    }else{
        queryglob = "*_{1,2}*{fastq,fastq.gz,fq,fq.gz}"
    }
    querydir = params.querydir
    procdir = params.outputdir
}else{
    //------ Default queryglob from local.config
    if(params.default_queryglob){
        default_queryglob = params.default_queryglob
    }else{
        default_queryglob = "*_{1,2}*{fastq,fastq.gz,fq,fq.gz}"
    }
    //------ querydir and queryglob
    if(params.querydir && params.queryglob){
        querydir = params.querydir
        queryglob = params.queryglob
    } else if(params.querydir && !params.queryglob) {
        querydir = params.querydir
        queryglob = "$default_queryglob"
    } else {
        switch(params.readtype) {
             // Case statement defined for 4 cases
             // Each case statement section has a break condition to exit the loop
             case "raw":
                querydir = params.rawdir
                if(params.queryglob){
                    queryglob = params.queryglob
                }else{
                    queryglob = "$default_queryglob"
                }
                break;
             case "fastp":
                querydir = params.procdir
                if(params.queryglob){
                    queryglob = params.queryglob
                }else{
                    queryglob = "fastp_$default_queryglob"
                }
                break;
             case "decont":
                querydir = params.procdir
                if(params.queryglob){
                    queryglob = params.queryglob
                }else{
                    queryglob = "decont_$default_queryglob"
                }
                break;
             default:
                //fastp
                querydir = params.procdir
                queryglob = "fastp_$default_queryglob"
                break;
          }
      //---------  outputdir
      if(params.outputdir == ""){
          outputdir = file(params.procdir)
      } else {
          outputdir = file(params.outputdir)
      }
  }
}
println ""
println "querydir : $querydir"
println "queryglob : $queryglob"
println "outputdir : $outputdir"
println ""

//============= Parse profilers ===========
Set profilers_expected = ['align', 'fastp', 'decont', 'kraken2', 'bracken', 'centrifuge', 'humann3', 'test', 'srst2', 'rgi']
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
include { ALIGN; MAPPED } from './modules/align'
include { FASTP } from './modules/fastp'
include { DECONT } from './modules/decontamination'
include { KRAKEN2; BRACKEN } from './modules/kraken2'
include { CENTRIFUGE } from './modules/centrifuge'
include { HUMANN3 } from './modules/humann3'
include { TEST } from './modules/test'
include { SRST2 } from './modules/srst2'
include { RGI } from './modules/rgi'

//=========== Kraken Preprocesses ===========
// https://groovy-lang.gitlab.io/101-scripts/basico/command_local-en.html
// For kraken, load database into ram.
// Remember to increase the space of /dev/shm
// Check /etc/fstab and add line
// none         /dev/shm            tmpfs       defaults,size=61G     0       0
// Not ideal for using in combination with --resume
if(params.krakenMMAP && profilers.contains('kraken2')){
    def result  = new StringBuilder()
    def error   = new StringBuilder()
    database = file(params.krakenDB)
    databaseDir = database.getParent()
    databaseName = database.getName()
    DB = "/dev/shm/database/${databaseName}"
    cmd0 = ["sh", "-c", "mkdir -p $DB"].execute()
    cmd0.consumeProcessOutput(result, error)
    println result
    println error
    println "Loading $params.krakenDB to $DB..."
    cmd = ["sh", "-c", "cp -r $params.krakenDB/*.k2d $DB"].execute()
    cmd.consumeProcessOutput(result, error)
    println result
    println error
    cmd.waitForOrKill(1000*60*10) //Wait 10 min before initiating the other steps
}

//---------- Main workflow ------------
workflow {
    //------ Get prefix ---------
    //// Closure after Channel.fromFilePairs() is used to parse the prefix
    //println "$querydir/**/$queryglob"
    querydirname = file(querydir).name
    ch_input = Channel.fromFilePairs("$querydir/**/$queryglob", flat: true, maxDepth: params.maxdepth) { file ->
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
    if(profilers.contains('align')){
        bwaIndex = file(params.bwaIndex)
        bwaIndexDir = bwaIndex.getParent()
        bwaIndexName = bwaIndex.getName()
        ALIGN(ch_reads, outputdir, params.bwaIndexDir, params.bwaIndexName)
        MAPPED(ALIGN.out.output.collect(), outputdir, params.bwaIndexName)
        //ALIGN.out.stdout.view { "ALIGN STDOUT:\n$it" }
        //MAPPED.out.stdout.view { "MAPPED STDOUT:\n$it" }
    }

    //------ Kraken + Bracken ---------
    if(profilers.contains('kraken2')){
        if(params.krakenMMAP){
            KRAKEN2(ch_reads, DB)
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
            ch_kraken = Channel.fromFilePairs("$params.procdir/**/kraken2/$dbname/*.{kraken2.report,kraken2.tax}",flat: true, size: 2, checkIfExists: true) { file -> file.getParent().getParent().getParent().name }
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
    if(profilers.contains('humann3')){
        HUMANN3(ch_reads, outputdir, params.humannDB_Uniref, params.humannDB_Chocophlan, params.humannDB_bt2Chocophlan)
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

//---------- Clean up ------------
// Removes database from ram to free up space.
workflow.onComplete {
    if(params.krakenMMAP && profilers.contains('kraken2')){
        println "Removing $DB..."
        cmd2 = ["sh", "-c", "rm -r $DB"].execute()
        cmd2.waitForOrKill(1000*10)
    }
}

//----------- TODO ------------
// Add customized log files
