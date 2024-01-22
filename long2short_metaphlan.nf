#!/usr/bin/env nextflow
// DSL 2 syntax
nextflow.enable.dsl=2

// help message
params.help = false
def helpMessage() {
    log.info"""
    ############################################################################
    Converting long reads to short reads, and perform metaphlan
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
Set profilers_expected = ['GEN_SR','BWA', 'bowtie', 'minimap', 'metaphlan2']
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
include { GEN_SR } from './modules/generate_SR'
include { BWA; BOWTIE; MINIMAP } from './modules/align'
include { METAPHLAN2 } from './modules/humann3'

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

    ch_reads = ch_preinput
    //ch_input = ch_preinput.filter { prefix, filenames ->   }

    //------- Display files ---------
    if(params.profilers.size() == 0 ) {
    	ch_reads.view()
    }

    //-------- Short Read Conversion --------
    if(profilers.contains('GEN_SR')){
        GEN_SR(ch_reads, outputdir)
        ch_sreads = GEN_SR.out.output
    }else{
        ch_sreads = ch_reads
    }
    //GEN_SR.out.stdout.view { "ALIGN STDOUT:\n$it" }
    //ch_sreads.view()

    //-------- Alignment ----------
    if(profilers.contains('BWA')){
        bwaIndex = file(params.bwaIndex)
        bwaIndexDir = bwaIndex.getParent()
        bwaIndexName = bwaIndex.getName()
        BWA(ch_sreads, outputdir, bwaIndexDir, bwaIndexName)
        ch_sam = BWA.out.sam
        //MAPPED(BWA.out.mapped_count.collect(), outputdir, bwaIndexName)
        //BWA.out.stdout.view { "ALIGN STDOUT:\n$it" }
        //MAPPED.out.stdout.view { "MAPPED STDOUT:\n$it" }
    }
    else if(profilers.contains('bowtie')){
        btIndex = file(params.btIndex)
        btIndexDir = btIndex.getParent()
        btIndexName = btIndex.getName()
        BOWTIE(ch_sreads, outputdir, btIndexDir, btIndexName)
        ch_sam = BOWTIE.out.sam
        //MAPPED(ALIGN.out.output.collect(), outputdir, bwaIndexName)
        //ALIGN.out.stdout.view { "ALIGN STDOUT:\n$it" }
        //MAPPED.out.stdout.view { "MAPPED STDOUT:\n$it" }
    }
    else if(profilers.contains('minimap')){
        // TBD: Add check database
        mmReadType = "sr"
        mmIndex = file(params.mmIndex)
        mmIndexDir = mmIndex.getParent()
        mmIndexName = mmIndex.getName()
        MINIMAP(ch_sreads, outputdir, mmReadType, mmIndexDir, mmIndexName)
        ch_sam = MINIMAP.out.sam
        MINIMAP.out.stdout.view { "MINIMAP STDOUT:\n$it" }
    }else{
        ch_sam = ch_sreads
    }

    //------ Metaphlan3 + Humann3 -------
    if(profilers.contains('metaphlan2')){
        METAPHLAN2(ch_sam, outputdir, params.metaphlan2DB_index_name, params.metaphlan2DB_pkl)
        METAPHLAN2.out.stdout.view()
    }
}
