`run.lrg.peaks` <-
function(k.biweight = 6, trans.method = "shiftedlog", add.par=10, root.dir = ".", 
        peak.dir, base.dir, lrg.dir, lrg.file = "lrg.peaks.RData",
        overwrite = FALSE, use.par.file = FALSE, par.file = "parameters.RData", 
        calc.all.peaks=TRUE){
    if(missing(base.dir)){base.dir <- paste(root.dir, "/Baseline_Corrected", sep="")}
    if(missing(peak.dir)){peak.dir <- paste(root.dir, "/All_Peaks", sep="")}
    if(missing(lrg.dir)){lrg.dir <- paste(root.dir, "/Large_Peaks", sep="")}
    if(use.par.file){
        load(paste(root.dir, "/", par.file, sep=""))
    }
    if(!file.exists(lrg.dir)){
        dir.create(lrg.dir)
    }
    if(!file.exists(paste(lrg.dir, "/", lrg.file, sep="")) || overwrite){
        file.list <- sub("_peaks\\.RData", "", list.files(peak.dir))
        for(i in list.files(peak.dir)){
            if(calc.all.peaks){
                load(paste(base.dir, "/", sub("_peaks.RData$", ".RData", i), sep=""))
                if(trans.method=="shiftedlog"){
                    peak.base$Amp <- log(peak.base$Amp+add.par-min(peak.base$Amp))
                }
                if(trans.method=="glog"){
                    peak.base$Amp <- log((peak.base$Amp+sqrt(add.par + peak.base$Amp^2))/2)
                }
                cent <- .biweight(peak.base$Amp, k.biweight)
                threshhold <- cent$center + k.biweight * cent$scale
                rm(peak.base)
            }
            load(paste(peak.dir, "/", i, sep=""))
            all.peaks$File <- factor(sub("_peaks\\.RData", "",i), file.list)
            if(i == list.files(peak.dir)[1]){
                if(calc.all.peaks){
                    lrg.peaks <- all.peaks[all.peaks$Max_hat >= threshhold,]
                } else {
                    lrg.peaks <- all.peaks
                }
            } else {
                if(calc.all.peaks){
                    lrg.peaks <- rbind(lrg.peaks, all.peaks[all.peaks$Max_hat >= threshhold,])
                } else {
                    lrg.peaks <- rbind(lrg.peaks, all.peaks)
                }
            }
            rm(all.peaks,threshhold)
        }
        save(lrg.peaks, file=paste(lrg.dir, "/", lrg.file, sep=""))
    } else {
        warning("File created by previous run of run.lrg.peaks exists and overwrite = FALSE; file not updated")
    }
}

