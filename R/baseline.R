`baseline` <-
function(spect, sm.par=1.1E-9, neg.pen=sqrt(pi/2), k.biweight=6, 
        max.iter=30, frac.changed=0.001, xvals=1:length(spect)){
    require(Matrix)
    L <- length(spect)
    s <- .biweight(spect, K=k.biweight)$scale/qnorm(0.75)
    neg.pen <- neg.pen/s  
    sm.fac <- L^4*sm.par/s

    D0 <- new("dgCMatrix", Dim=as.integer(c(L,L)), p=as.integer(cumsum(c(0,3,4,rep(5,L-4),4,3))), 
        i=as.integer(c(0:2,0:3,as.vector(outer(0:4,0:(L-5),"+")), (L-4):(L-1), (L-3):(L-1))), 
        x=c(1,-2,1,-2,5,-4,1,rep(c(1,-4,6,-4,1),L-4),1,-4,5,-2,1,-2,1))
    
    bd <- rep(median(spect),L)
    bd0 <- spect

    indicator <- (bd>spect)
    indicator0 <- rep(TRUE,L)
    changed <- c()
    iter <- 0
    while(sum(xor(indicator,indicator0)) > frac.changed * length(spect) && iter < max.iter){
        indicator0 <- indicator
        M <- Matrix(1/2/sm.fac + neg.pen*spect*indicator/sm.fac)
        D1 <- D0 + new("dgCMatrix", Dim=as.integer(c(L,L)), x=neg.pen*indicator/sm.fac, 
            p=as.integer(0:L), i=as.integer(0:(L-1)))
  
        bd <- as.numeric(solve(D1,M))
        indicator <- (bd > spect)
        changed <- c(changed,sum(xor(indicator,indicator0)))
        iter <- iter + 1
    }
    if(sum(xor(indicator,indicator0)) > frac.changed * length(spect)){
        warning("Iteration limit reached without convergence to specified tolerance.")
    }
    list(baseline = bd, noise = s, iter = iter, changed = changed)
}

