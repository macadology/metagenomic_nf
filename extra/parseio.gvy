//def parseio(params) {
if(params.general){
    //Checks if local.config is a placeholder or specific to project.
    if(params.outputdir == "" || params.querydir == ""){
        throw new Exception("Error: You are running in general mode. Please specify a query (--querydir) and output (--outputdir) directory")
    }
    //------ queryglob, querydir, outputdir
    if(params.queryglob){
        queryglob = params.queryglob
    }else{
        queryglob = "*_{1,2}*{fastq,fastq.gz,fq,fq.gz}"
    }
    querydir = params.querydir
    outputdir = params.outputdir
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
    }
    //---------  outputdir
    if(!params.outputdir){
        println params.procdir.getClass()
        outputdir = params.procdir
    } else {
        outputdir = params.outputdir
    }
}
//    return [querydir, queryglob, outputdir]
//}
