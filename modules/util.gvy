def parse_general(params_outputdir, params_querydir, params_queryglob, default_queryglob) {
    //---
    // Checks directories and query glob if the general = True
    //---
    //Checks if local.config is a placeholder or specific to project.
    if(outputdir == "" || querydir == ""){
        log.error "Error: You are running in general mode. Please specify a query (--querydir) and output (--outputdir) directory"
        exit 0
    }
    //------ queryglob, querydir, outputdir
    if(params_queryglob){
        queryglob = params_queryglob
    }else{
        queryglob = default_queryglob
    }
    querydir = params_querydir
    procdir = params_outputdir
}

def parse_defaults(params_default_queryglob, params_queryglob, TBA) {
    //------ Default queryglob from local.config
    if(params.default_queryglob){
        default_queryglob = params_default_queryglob
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
