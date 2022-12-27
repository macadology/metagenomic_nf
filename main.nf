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
                                show matched reads. Options -> 'BWA', 'fastp', 'decont', 'kraken2', 'bracken', 'centrifuge', 'humann3', 'srst2', 'rgi'
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
        nextflow run main.nf -profile batch -plugins nf-amazon --bucket-dir s3://jon-nextflow-work --profilers humann3

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

//=========== Load functions ============
// GroovyShell shell = new GroovyShell()
// def funcs = shell.parse(new File('./extra/functions.gvy'))
// (querydir, queryglob, outputdir) = funcs.parseio(params)

evaluate(new File("./extra/parseio.gvy"))

println ""
println "querydir : $querydir"
println "queryglob : $queryglob"
println "query : $querydir/**/$queryglob"
println "outputdir : $outputdir"
println "database: $params.database"
println ""

//============= Parse profilers ===========
Set profilers_expected = ['BWA', 'bowtie', 'minimap', 'sam2fastq', 'fastp', 'decont', 'kraken2', 'bracken', 'centrifuge', 'humann3', 'metaphlan', 'test', 'srst2', 'rgi']
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
include { BWA; BOWTIE; MINIMAP; SAM2FASTQ; MAPPED } from './modules/align'
include { FASTP } from './modules/fastp'
include { DECONT } from './modules/decontamination'
include { KRAKEN2; BRACKEN } from './modules/kraken2'
include { CENTRIFUGE } from './modules/centrifuge'
include { HUMANN3; METAPHLAN } from './modules/humann3'
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
    //------ Get prefix and files ---------
    querydirname = new File(querydir).name
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
        //ch_reads.view()
        println "testDatabase: $params.testDatabase"
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
