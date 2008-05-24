`run.peaks` <-
function(add.par=10, trans.method="shiftedlog", root.dir=".", base.dir, 
        peak.dir, overwrite=FALSE, use.par.file=FALSE, par.file = "parameters.RData", 
        num.pts=5, R2.thresh=0.98, oneside.min=1, peak.method="parabola", 
        calc.all.peaks=TRUE, numsds=4){
    fail <- 0
    if(missing(base.dir)){base.dir <- paste(root.dir, "/Baseline_Corrected", sep="")}
    if(missing(peak.dir)){peak.dir <- paste(root.dir, "/All_Peaks", sep="")}
    if(use.par.file){
        load(paste(root.dir, "/", par.file, sep=""))
    }
    if(!file.exists(peak.dir)){
        dir.create(peak.dir)
    }
    for(j in list.files(base.dir)){
        if(!file.exists(paste(peak.dir, "/", sub("\\.RData$", "_peaks.RData", j), sep="")) || 
                overwrite){
            load(paste(base.dir, "/", j, sep=""))
            names(peak.base) <- c("Freq","LogAmp")
            if(trans.method=="shiftedlog"){
                peak.base$LogAmp <- log(peak.base$LogAmp+add.par-min(peak.base$LogAmp))
            }
            if(trans.method=="glog"){
                peak.base$LogAmp <- log((peak.base$LogAmp+sqrt(add.par + peak.base$LogAmp^2))/2)
            }
            thresh <- ifelse(calc.all.peaks, -Inf, mean(peak.base$LogAmp) + numsds * sd(peak.base$LogAmp))
            peak.base <- peak.base[order(peak.base$Freq),]
            numsplit <- (dim(peak.base)[1] %/% 100000) + 1
            splits <- (1:(numsplit-1)) * floor(dim(peak.base)[1]/numsplit)
            splits <- apply(cbind(c(1,splits-(num.pts-oneside.min-1)),
                c(splits+(num.pts-oneside.min-1),dim(peak.base)[1])), 1, list)
            splits <- lapply(splits, function(x){seq(x[[1]][1],x[[1]][2])})
            for(i in 1:numsplit){
                all.peaks.tmp <- locate.peaks(peak.base[splits[[i]],], num.pts, 
                    R2.thresh, oneside.min, peak.method, thresh)
                save(all.peaks.tmp, file=paste(peak.dir, "/", sub("\\.RData", 
                    paste("_peaks", i, ".RData", sep=""), j), sep=""))
                rm(all.peaks.tmp)
            }
            for(i in numsplit:1){
                load(paste(peak.dir, "/", sub("\\.RData", paste("_peaks", i, ".RData", 
                    sep=""), j), sep=""))
                if(i == numsplit){
                    all.peaks <- all.peaks.tmp
                } else {
                    all.peaks <- rbind(all.peaks.tmp,all.peaks)
                }
                rm(all.peaks.tmp)
                file.remove(paste(peak.dir, "/", sub("\\.RData", paste("_peaks", i, 
                    ".RData", sep=""), j), sep=""))
            }
            all.peaks <- unique(all.peaks)
            rownames(all.peaks) <- 1:dim(all.peaks)[1]
            save(all.peaks, file=paste(peak.dir, "/", sub("\\.RData$", 
                "_peaks.RData", j), sep=""))
        } else {
            fail <- fail + 1
        }
    }
    if(fail) {
        warning(paste(fail, "peak file(s) already existed and overwrite = FALSE; those file(s) not updated"))
    }
}

