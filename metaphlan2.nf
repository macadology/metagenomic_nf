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
Set profilers_expected = ['BWA', 'bowtie', 'minimap', 'metaphlan2']
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
include { BWA; BOWTIE; MINIMAP} from './modules/align'
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

    //-------- Alignment ----------
    if(profilers.contains('BWA')){
        bwaIndex = file(params.bwaIndex)
        bwaIndexDir = bwaIndex.getParent()
        bwaIndexName = bwaIndex.getName()
        BWA(ch_reads, outputdir, bwaIndexDir, bwaIndexName)
        ch_sam = BWA.out.sam
        //MAPPED(BWA.out.mapped_count.collect(), outputdir, bwaIndexName)
        //BWA.out.stdout.view { "ALIGN STDOUT:\n$it" }
        //MAPPED.out.stdout.view { "MAPPED STDOUT:\n$it" }
    }else{
        ch_sam = ch_reads
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
        ch_sam = ch_reads
    }

    if(profilers.contains('minimap')){
        // TBD: Add check database
        mmReadType = "sr"
        mmIndex = file(params.mmIndex)
        mmIndexDir = mmIndex.getParent()
        mmIndexName = mmIndex.getName()
        MINIMAP(ch_reads, outputdir, mmReadType, mmIndexDir, mmIndexName)
        ch_sam = MINIMAP.out.sam
        MINIMAP.out.stdout.view { "MINIMAP STDOUT:\n$it" }
    }else{
        ch_sam = ch_reads
    }

    //------ Metaphlan3 + Humann3 -------
    if(profilers.contains('metaphlan2')){
        METAPHLAN2(ch_sam, outputdir, params.metaphlan2DB_index_name, params.metaphlan2DB_pkl)
        METAPHLAN2.out.stdout.view()
    }
}
