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
Set profilers_expected = ['fastp', 'decont', 'bowtie', 'metaphlan', 'humann3']
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
include { HUMANN3 } from '../modules/humann3'


process STAGEIN {
    maxForks 3
    stageInMode 'copy'
    executor 'local'

    input:
    tuple val(prefix), path(reads)

    output:
    tuple val(prefix), path(reads), emit: output

    script:
    """
    """
}

process CLEANFILES {
    cache 'lenient'
    cpus 4
    memory 16.GB
    queue 'express'
    time '1hour'

    input:
    val(files_input)

    output:
    val(1), emit: IS_CLEAN
    stdout emit: stdout

    script:
    """
    for file in ${files_input}; do
      # Remove cruff added by Nextflow
      if [ -e \$file ]; then
        # Log some info about the file for debugging purposes
        echo "cleaning \$file"
        # stat \$file
        # Get file info: size, access and modify times
        size=`stat --printf="%s" \$file`
        atime=`stat --printf="%X" \$file`
        mtime=`stat --printf="%Y" \$file`

        # Make the file size 0 and set as a sparse file
        > \$file
        truncate -s \$size \$file
        # Reset the timestamps on the file
        touch -a -d @\$atime \$file
        touch -m -d @\$mtime \$file
      fi
    done
    """
}

process TEST {
    // Copy files from DATA to scratch
    cpus 4
    memory 16.GB
    queue 'express'
    time '1hour'
    maxForks 3
    scratch "$HOME/scratch/temp"
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }

    input:
    tuple val(prefix), path(reads)

    output:
    publishDir "/home/users/astar/gis/teojyj/DATA/test", mode: 'move'
    tuple val(prefix), path("${prefix}.txt"), emit: output
    stdout emit: stdout

    script:
    //String indexname = "${bwaIndex.baseName}_${bwaIndex.extension}"
    """
    echo "Executing test on ${prefix}"
    cat ${reads} > ${prefix}.txt
    echo "haha" > ${prefix}_haha.txt
    sleep 10
    """
}

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

    //https://pirl.unc.edu/blog/tricking-nextflows-caching-system-to-drastically-reduce-storage-usage
    if(profilers.contains('humann3')){
        STAGEIN(ch_input)
        ch_staged = STAGEIN.out
        HUMANN3(ch_staged, outputdir, params.humannDB_index, params.humannDB_Uniref, params.humannDB_Chocophlan, params.humannDB_bt2Chocophlan)
        HUMANN3.out.stdout.view()
        //TEST(ch_staged)
        STAGEIN.output
            .join(HUMANN3.out.output, by: [0])
            .flatten()
            .filter{ it =~ /.*fq.gz/ }
            .set{ ch_done }
        CLEANFILES(ch_done)
        CLEANFILES.out.stdout.view()
        //if( params.delete_intermediates ) {
        //    clean_work_files(ch_done)
        //}
    }
}
