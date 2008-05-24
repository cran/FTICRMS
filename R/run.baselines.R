`run.baselines` <-
function(root.dir = ".", raw.dir, base.dir, overwrite=FALSE, use.par.file = FALSE, 
        par.file = "parameters.RData", sm.fac=10^15, neg.pen=sqrt(pi/2), max.iter=30, 
        frac.changed=0.001, binsize=64){
    fail <- 0
    if(missing(base.dir)){base.dir <- paste(root.dir, "/Baseline_Corrected", sep="")}
    if(missing(raw.dir)){raw.dir <- paste(root.dir, "/Raw_Data", sep="")}
    if(use.par.file){
        load(paste(root.dir, "/", par.file, sep=""))
    }
    if(!file.exists(base.dir)){
        dir.create(base.dir)
    }
    for(i in list.files(raw.dir)){
        if(!file.exists(paste(base.dir, "/", sub("\\.txt$", ".RData", i), sep="")) || 
                overwrite){
            if(regexpr(",", readLines(paste(raw.dir, "/", i, sep=""), n=1)) != -1){ 
                peak.base <- read.csv(paste(raw.dir, "/", i, sep=""), header=FALSE)
            } else {
                peak.base <- read.table(paste(raw.dir, "/", i, sep=""), header=FALSE)
            }
            peak.base[,2] <- peak.base[,2] - baseline(peak.base[,2], sm.fac, 
                neg.pen, max.iter, frac.changed, binsize)[[1]]
            names(peak.base) <- c("Freq", "Amp")
            save(peak.base, file=sub("\\.txt$",".RData", paste(base.dir, "/", i, sep="")))
            rm(peak.base)
        } else {
            fail <- fail + 1
        }
    }
    if(fail) {
        warning(paste(fail, "baseline file(s) already existed and overwrite = FALSE; those file(s) not updated"))
    }
}

