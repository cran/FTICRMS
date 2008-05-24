`run.cluster.matrix` <-
function(pre.align = FALSE, align.method = "spline", trans.method="shiftedlog", 
        add.par = 10, lrg.only = TRUE, masses = NULL, isotope.dist = 7, 
        root.dir = ".", base.dir, peak.dir, lrg.dir, lrg.file = "lrg.peaks.RData", 
        overwrite = FALSE, use.par.file = FALSE, par.file = "parameters.RData"){
    if(missing(base.dir)){base.dir <- paste(root.dir, "/Baseline_Corrected", sep="")}
    if(missing(peak.dir)){peak.dir <- paste(root.dir, "/All_Peaks", sep="")}
    if(missing(lrg.dir)){lrg.dir <- paste(root.dir, "/Large_Peaks", sep="")}
    if(use.par.file){
        load(paste(root.dir, "/", par.file, sep=""))
    }
    if(is.null(masses) && !lrg.only){
        stop("Either lrg.only must be TRUE or masses must be defined (or both)")
    }
    load(paste(lrg.dir, "/", lrg.file, sep=""))
    if(!exists("clust.mat") || overwrite){
        if(is.na(amps) && align.method!="none"){
            align.method <- "none"
            warning("Parameter 'align.method' changed to 'none' because no peaks appeared in all samples")
            if(use.par.file){
                tmp <- extract.pars(root.dir, par.file)
                tmp$align.method <- "none"
                do.call(make.par.file, tmp)
                load(paste(root.dir, "/", par.file, sep=""))
            }
        }
        if(!identical(pre.align, FALSE)){
            if(!identical(class(pre.align),"list")){
                pre.align <- list(targets=0, actual=data.frame(-pre.align))
                if(use.par.file){
                    tmp <- extract.pars(root.dir, par.file)
                    tmp$pre.align <- pre.align
                    do.call(make.par.file, tmp)
                    load(paste(root.dir, "/", par.file, sep=""))
                }
            }
            pre.spl <- do.call(.alignment, pre.align)
        }
        if(align.method=="spline"){
            require(splines)
            targets <- sapply(centers,mean)
            act.list <- lapply(1:dim(centers)[1], function(x){as.numeric(centers[x,])})
            spl <- lapply(act.list, function(x){interpSpline(x,targets)})
            names(spl) <- rownames(centers)
            for(i in names(spl)){
                lrg.peaks$Center_hat[lrg.peaks$File==i] <-
                    predict(spl[[i]],lrg.peaks$Center_hat[lrg.peaks$File==i])$y
            }
        }

        peaks.clustered <- .cluster.peaks(lrg.peaks)
        peaks.clustered <- .break.dups(peaks.clustered)
        clust.mat <- .cluster.matrix(peaks.clustered)
        rm(peaks.clustered)
        if(!is.null(masses)){
            wh <- .get.sp.masses(rownames(clust.mat), masses, isotope.dist)
            clust.mat <- clust.mat[wh,]
        }
        if(!lrg.only){
            add.names <- .get.sp.masses(rownames(clust.mat), masses, isotope.dist, inds=FALSE)
            add.rows <- matrix(NA, ncol=dim(clust.mat)[2], nrow=length(add.names))
            add.rows <- data.frame(add.rows)
            rownames(add.rows) <- add.names
            colnames(add.rows) <- colnames(clust.mat)
            clust.mat <- rbind(clust.mat, add.rows)
            clust.mat <- clust.mat[order(sapply(rownames(clust.mat), function(x)
                mean(as.numeric(strsplit(x, ",")[[1]])))),]
            tmprows <- do.call(rbind, strsplit(rownames(clust.mat),","))
            wh <- tmprows[1:(dim(tmprows)[1]-1),2] < tmprows[2:dim(tmprows)[1],1]
            clust.mat <- clust.mat[(c(TRUE,wh) & c(wh,TRUE)) | apply(!is.na(clust.mat),1,any),]
            rm(add.rows,add.names,tmprows,wh)
        }
        num.sig <- apply(!is.na(clust.mat), 1, sum)
        ran <- matrix(sapply(strsplit(rownames(clust.mat),","),as.numeric),
            ncol=2,byrow=TRUE)
        for(f in names(clust.mat)){
            load(paste(peak.dir, "/", f, "_peaks.RData", sep=""))
            if(!identical(pre.align, FALSE)){
                if(identical(class(pre.spl[[f]]),"lm")){
                    all.peaks$Center_hat <- predict(pre.spl[[f]], data.frame(X=all.peaks$Center_hat))
                } else {
                    all.peaks$Center_hat <- predict(pre.spl[[f]], all.peaks$Center_hat)$y
                }
            }
            if(align.method=="spline"){
                all.peaks$Center_hat <- predict(spl[[f]], all.peaks$Center_hat)$y
            }
            locs <- which(is.na(clust.mat[,f]))
            endpts <- as.numeric(t(ran[locs,,drop=FALSE]))
            ints <- findInterval(all.peaks$Center_hat, endpts, rightmost.closed=TRUE)
            ints <- (ints * (ints %% 2)+1)/2
            all.peaks <- data.frame(Int = ints[ints>.75], Max_hat=all.peaks$Max_hat[ints>.75])
            if(dim(all.peaks)[1]){
                all.peaks <- unlist(by(all.peaks$Max_hat,all.peaks$Int,max))
                clust.mat[locs[as.numeric(names(all.peaks))],f] <- all.peaks
            }
            rm(all.peaks)

            if(any(is.na(clust.mat[,f]))){
                load(paste(base.dir, "/", f, ".RData", sep=""))
                names(peak.base) <- c("Freq","LogAmp")
                if(trans.method == "shiftedlog"){
                    peak.base$LogAmp <- log(peak.base$LogAmp+add.par-min(peak.base$LogAmp))
                }
                if(trans.method=="glog"){
                    peak.base$LogAmp <- log((peak.base$LogAmp+sqrt(add.par + peak.base$LogAmp^2))/2)
                }
                if(!identical(pre.align, FALSE)){
                    if(identical(class(pre.spl[[f]]),"lm")){
                        peak.base$Freq <- predict(pre.spl[[f]], data.frame(X=peak.base$Freq))
                    } else {
                        peak.base$Freq <- predict(pre.spl[[f]], peak.base$Freq)$y
                    }
                }
                if(align.method=="spline"){
                    peak.base$Freq <- predict(spl[[f]], peak.base$Freq)$y
                }
                peak.base <- peak.base[order(peak.base$Freq),]
                locs <- which(is.na(clust.mat[,f]))
                endpts <- as.numeric(t(ran[locs,,drop=FALSE]))
                ints <- findInterval(peak.base$Freq, endpts, rightmost.closed=TRUE)
                ints <- (ints * (ints %% 2)+1)/2
                all.peaks <- data.frame(Int = ints[ints>.75], Max_hat=peak.base$LogAmp[ints>.75])
                if(dim(all.peaks)[1]){
                    all.peaks <- unlist(by(all.peaks$Max_hat,all.peaks$Int,max))
                    clust.mat[locs[as.numeric(names(all.peaks))],f] <- all.peaks
                }
                if(any(is.na(clust.mat[,f]))){
                    locs <- which(is.na(clust.mat[,f]))
                    ends <- findInterval(apply(ran[locs,,drop=FALSE],1,mean), peak.base$Freq,
                        all.inside=TRUE)
                    ends2 <- apply(abs(cbind(peak.base$Freq[ends],peak.base$Freq[ends+1]) -
                        apply(ran[locs,,drop=FALSE],1,mean)), 1, which.min)
                    clust.mat[locs,f] <- peak.base$LogAmp[ends+ends2-1]

                }
                rm(peak.base)
            }
        }
        save(lrg.peaks,clust.mat,amps,centers,num.sig, file=paste(lrg.dir, "/", 
            lrg.file, sep=""))
    } else {
        warning("File created by previous run of run.cluster.matrix exists and overwrite = FALSE; file not updated")
    }
}
