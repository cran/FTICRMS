`display.tests` <-
function (sig.rows = "all", summ = anova, tests, form = parameter.list$form){
    if (missing(tests)) {
        if (identical(sig.rows, "all")) {
            sig.rows <- 1:dim(sigs)[1]
        }
        tests <- match(rownames(sigs)[sig.rows], rownames(clust.mat))
        if (any(is.na(tests))) {
            tests <- tests[!is.na(tests)]
            warning("Nonexistent row requested from clust.mat")
        }
    }
    ret <- lapply(tests, function(x) {
        dat <- data.frame(Y=t(clust.mat[x,,drop=FALSE]), unique(parameter.list$covariates))
        colnames(dat)[1] <- "Y"
        lm(form, dat)
    })
    if(!identical(summ,"none")){
        ret <- lapply(ret, summ)
    }
    names(ret) <- rownames(clust.mat)[tests]
    ret
}
