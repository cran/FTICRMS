`extract.pars` <-
function(root.dir=".", par.file = "parameters.RData"){
    load(paste(root.dir, "/", par.file, sep=""))
    ret <- list()
    for(i in setdiff(ls(),c("i","ret"))){
        ret[[i]] <- get(i)
    }
    ret
}

