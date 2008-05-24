`baseline` <-
function(spect, sm.fac=10^15, neg.pen=sqrt(pi/2), 
        max.iter=30, frac.changed=0.001, binsize=64, xvals=1:length(spect)){
    require(Matrix)
    L <- length(spect)
    P <- L %/% binsize
    s <- quantile(apply(matrix(spect[1:(P*binsize)],nrow=binsize),2,sd),0.25)
    neg.pen <- neg.pen/s  
    sm.fac <- sm.fac/s

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

