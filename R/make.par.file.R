`make.par.file` <-
function(covariates, form, root.dir=".", par.file="parameters.RData", ...){
    default.vals <- list(
        add.norm = TRUE,
        add.par = 10,
        align.method = "spline",
        base.dir = paste(root.dir, "/Baseline_Corrected", sep=""),
        calc.all.peaks = TRUE,
        clust.constant = 10,
        clust.method = "ppm",
        cor.thresh = 0.8,
        FDR = 0.1,
        frac.changed = 0.001,
        isotope.dist = 7,
        k.biweight = 6,
        lrg.dir = paste(root.dir, "/Large_Peaks", sep=""), 
        lrg.file = "lrg.peaks.RData",
        lrg.only = TRUE,
        masses = NULL,
        max.iter = 30,
        neg.pen = sqrt(pi/2),
        normalization = "common",
        num.pts = 5,
        oneside.min = 1,
        overwrite = FALSE,
        par.file = "parameters.RData",
        peak.dir = paste(root.dir, "/All_Peaks", sep=""),
        peak.method = "parabola",
        pre.align = FALSE,
        pval.fcn = "default",
        R2.thresh = 0.98,
        raw.dir = paste(root.dir, "/Raw_Data", sep=""),
        repl.method = max,
        res.dir = paste(root.dir, "/Results", sep=""),
        res.file = "analyzed.RData",
        root.dir = ".",
        sm.par = 1.1E-9,
        trans.method = "shiftedlog",
        use.t.test = TRUE
    )
    new.vals <- c(list(...), root.dir=root.dir, par.file=par.file)
    for(i in names(default.vals)){
        if(i %in% names(new.vals)){
            assign(i, new.vals[[i]])
        } else {
            assign(i, default.vals[[i]])        
        }
    }
    tmp <- !(names(new.vals) %in% names(default.vals))
    if(sum(tmp)==1){
        warning(paste(names(new.vals)[tmp], "is not a valid parameter name"))
    } else if(sum(tmp)==2){
        warning(paste(paste(names(new.vals)[tmp], collapse=" and "), 
            "are not valid parameter names"))
    } else if(sum(tmp)>2){
        names(new.vals)[tmp][sum(tmp)] <- paste("and", names(new.vals)[tmp][sum(tmp)])
        warning(paste(paste(names(new.vals)[tmp], collapse=", "), 
            "are not valid parameter names"))
    }
    save(list=setdiff(ls(), c("default.vals", "new.vals", "tmp", "i")), 
        file=paste(root.dir, "/", par.file, sep=""))
}

