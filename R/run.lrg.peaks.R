`run.lrg.peaks` <-
function(trans.method = "shiftedlog", add.par = 0, subtract.base = FALSE, root.dir = ".",
        peak.dir, base.dir, lrg.dir, lrg.file = "lrg_peaks.RData",
        overwrite = FALSE, use.par.file = FALSE, par.file = "parameters.RData",
        calc.all.peaks = FALSE, gengamma.quantiles = TRUE, peak.thresh = 3.798194){
    if(missing(base.dir)){base.dir <- paste(root.dir, "/Baselines", sep="")}
    if(missing(peak.dir)){peak.dir <- paste(root.dir, "/All_Peaks", sep="")}
    if(missing(lrg.dir)){lrg.dir <- paste(root.dir, "/Large_Peaks", sep="")}
    if(use.par.file){
        load(paste(root.dir, "/", par.file, sep=""))
    }
    if(!file.exists(lrg.dir)){
        dir.create(lrg.dir)
    }
    zeros <- 0
    if(!file.exists(paste(lrg.dir, "/", lrg.file, sep="")) || overwrite){
        file.list <- sub("_peaks\\.RData", "", list.files(peak.dir))
        for(i in list.files(peak.dir)){
            if(calc.all.peaks){
                load(paste(base.dir, "/", i, ".RData", sep=""))
                names(spect) <- c("Freq","LogAmp")
                if(subtract.base){
                    spect$LogAmp <- spect$LogAmp-spect.base
                } else if(any(spect$LogAmp == 0)){
                    zeros <- zeros + 1
                    spect$LogAmp[spect$LogAmp == 0] <- min(spect$LogAmp[spect$LogAmp > 0])
                }
                if(trans.method=="shiftedlog"){
                    if(subtract.base){
                        add.par <- add.par - min(spect$LogAmp)
                    }
                    spect$LogAmp <- log(spect$LogAmp+add.par)
                }
                if(trans.method=="glog"){
                    spect$LogAmp <- log((spect$LogAmp+sqrt(add.par + spect$LogAmp^2))/2)
                }
                if(gengamma.quantiles){
                    thresh <- log(peak.thresh*spect.base)
                } else {
                    thresh <- .biweight.FTICRMS(spect$LogAmp, K=3*peak.thresh/2)
                    thresh <- thresh$center + 3*peak.thresh/2*thresh$scale
                }
                rm(spect,spect.base)
            }
            load(paste(peak.dir, "/", i, sep=""))
            all.peaks$File <- factor(sub("_peaks\\.RData", "",i), file.list)
            if(i == list.files(peak.dir)[1]){
                if(calc.all.peaks){
                    thresh <- approx(spect$Freq, thresh, all.peaks$Center_Hat)$y
                    lrg.peaks <- all.peaks[all.peaks$Max_hat >= thresh,]
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
            rm(all.peaks)
            if(calc.all.peaks){
                rm(thresh)
            }
        }
        if(zeros){
            warning(paste(zeros, ifelse(zeros==1, "spectrum", "spectra"), "had one or more zero entries for amplitude; those entries replaced by minimum value of amplitude"))
        }
        save(lrg.peaks, file=paste(lrg.dir, "/", lrg.file, sep=""))
    } else {
        warning("File created by previous run of run.lrg.peaks exists and overwrite = FALSE; file not updated")
    }
}

