`run.analysis` <-
function(form, covariates, FDR=0.1, normalization="common", add.norm=TRUE, 
        repl.method=max, use.t.test = TRUE, lrg.only = TRUE, masses = NULL, 
        isotope.dist = 7, root.dir=".", lrg.dir, lrg.file = "lrg.peaks.RData", 
        res.dir, res.file = "analyzed.RData", overwrite = FALSE, 
        use.par.file = FALSE, par.file = "parameters.RData", ...){
    if(missing(res.dir)){res.dir <- paste(root.dir, "/Results", sep="")}
    if(missing(lrg.dir)){lrg.dir <- paste(root.dir, "/Large_Peaks", sep="")}
    if(use.par.file){
        load(paste(root.dir, "/", par.file, sep=""))
        parameter.list <- extract.pars(root.dir, par.file)
    } else {
        parameter.list <- NA
    }
    if(is.null(masses) && !lrg.only){
        stop("Either lrg.only must be TRUE or masses must be defined (or both)")
    }
    if(!file.exists(res.dir)){
        dir.create(res.dir)
    }
    if(!file.exists(paste(res.dir, "/", res.file, sep="")) || overwrite){
        load(paste(lrg.dir, "/", lrg.file, sep=""))
        norm <- switch(normalization,
            none = dim(clust.mat)[2],
            postbase = by(lrg.peaks$Max_hat,lrg.peaks$File,mean),
            postrepl = dim(clust.mat)[2],
            common = apply(amps,1,mean)
        )
        if(add.norm){
            clust.mat <- clust.mat - rep(norm, each=dim(clust.mat)[1]) + mean(norm)
        } else {
            clust.mat <- as.matrix(clust.mat) %*% diag(mean(1/norm)/norm)
        }
        if(!identical(repl.method,"none")){
            if(is.character(repl.method)){
                repl.method <- get(repl.method)
            }
            clust.mat <- do.call(cbind, by(as.data.frame(t(clust.mat)), covariates$subj, function(x)apply(x,2,repl.method)))
            if(length(cols <- grep("^[^:]*$",terms(form)@term.labels,value=TRUE))==1){
                covariates <- data.frame(c(by(covariates[,cols], covariates$subj, unique)))
                colnames(covariates) <- cols
            } else {
                covariates <- do.call(rbind, by(covariates[,cols], 
                    covariates$subj, unique))
            }
        }
        if(normalization == "postrepl"){
            norm <- apply(clust.mat, 2, mean)
            if(add.norm){
                clust.mat <- clust.mat - rep(norm, each=dim(clust.mat)[1]) + mean(norm)
            } else {
                clust.mat <- as.matrix(clust.mat) %*% diag(mean(1/norm)/norm)
            }
        }
        
        if(!is.null(masses)){
            wh <- .get.sp.masses(names(num.sig), masses, isotope.dist)
            clust.mat <- clust.mat[wh,]
            num.sig <- num.sig[wh]
        }
    
        p.value <- rep(0, dim(clust.mat)[1])
        form <- update(form, Y~.)
        if(use.t.test){
            Delta = rep(0, dim(clust.mat)[1])
            for(i in 1:length(p.value)){
                tmpdat <- data.frame(Y=t(clust.mat[i,,drop=FALSE]), covariates)
                colnames(tmpdat)[1] <- "Y"
                tmp <- t.test(form, dat=tmpdat, ...)
                p.value[i] <- tmp$p.value
                Delta[i] <- diff(tmp$estimate)
            }
            which.sig <- data.frame(Delta, p.value, num.sig, ord=1:length(p.value))
        } else {
            for(i in 1:length(p.value)){
                tmpdat <- data.frame(Y=t(clust.mat[i,,drop=FALSE]), covariates)
                colnames(tmpdat)[1] <- "Y"
                tmp <- summary(lm(form, dat=tmpdat, ...))$fstatistic
                p.value[i] <- pf(tmp[1],tmp[2],tmp[3],lower.tail=FALSE)
            }
            which.sig <- data.frame(Delta = NA, p.value, num.sig, ord=1:length(p.value))
        }
        which.sig <- which.sig[order(which.sig$p.val),]
        tmp <- matrix(NA, ncol=max(num.sig), nrow=dim(which.sig)[1])
        colnames(tmp) <- paste("S", 1:dim(tmp)[2], sep="")
        for(k in 1:dim(tmp)[2]){
            inds <- which(which.sig$num.sig >= k)
            tmp[which.sig$ord[inds],k] <- 0
            sp <- .benj.hoch(which.sig$p.val[inds], FDR)
            if(sp > 0){
                tmp[which.sig$ord[inds][1:sp],k] <- 1
            }
            rm(sp,inds)
        }
        if(!lrg.only){
            tmp <- cbind(0, tmp)
            colnames(tmp)[1] <- "S0"
            sp <- .benj.hoch(which.sig$p.val, FDR)
            if(sp > 0){
                tmp[which.sig$ord[1:sp],1] <- 1
            }
        }

        which.sig <- which.sig[order(which.sig$ord),c("Delta","p.value","num.sig")]
        which.sig <- data.frame(which.sig, tmp)
        sigs <- which.sig[apply(which.sig[,-(1:3),drop=FALSE],1,any, na.rm=TRUE),]
        sigs <- sigs[,apply(sigs,2,any,na.rm=TRUE)]
        min.FDR <- sapply(1:max(num.sig), function(x){
            tmp <- sort(which.sig$p.value[which.sig$num.sig>=x])
            min(length(tmp)*tmp/(1:length(tmp)))
        })
        names(min.FDR) <- 1:length(min.FDR)

        save(amps,centers,clust.mat,min.FDR,sigs,which.sig,parameter.list, 
            file=paste(res.dir, "/", res.file, sep=""))
    } else {
        warning("Results file exists and overwrite = FALSE; no results file created")
    }
}
