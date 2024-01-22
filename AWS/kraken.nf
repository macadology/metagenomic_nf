#!/usr/bin/env nextflow
// DSL 2 syntax
nextflow.enable.dsl=2

// help message
params.help = false
def helpMessage() {
    log.info"""
    Running Kraken on AWS

    Example code:
    nextflow kraken.nf --querydir [folder] --queryglob "chopper*.fastq.gz" --size 1 --outputdir [folder] --krakenDB /home/ec2-user/mnt/databases/k2_standard_08gb_20231009 --krakenReadlength 100 --krakenMMAP --krakenKeepOutput true --krakenThreads 12 --profilers kraken2,bracken

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
Set profilers_expected = ['kraken2']
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
include { KRAKEN2; BRACKEN } from './modules/kraken2'

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
    ch_preinput = Channel.fromFilePairs("$querydir/**/$queryglob", flat: false, size:params.size, maxDepth: params.maxdepth) { file ->
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

    ch_input = ch_preinput
    //ch_input = ch_preinput.filter { prefix, filenames ->   }

    //------- Display files ---------
    if(params.profilers.size() == 0 ) {
    	ch_input.view()
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
            if(params.outputdir == "" || params.querydir == ""){
                throw new Exception("Error: You are running bracken only. Please specify a query (--querydir <dir>), queryglob that include the database name (--queryglob <dbname>/*.kraken2.report), and output (--outputdir) directory")
            }
            ch_kraken = ch_input
        }
        BRACKEN(ch_kraken, outputdir, params.krakenDB)
        //BRACKEN.out.stdout.view { "BRACKEN STDOUT:\n$it" }
    }
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
