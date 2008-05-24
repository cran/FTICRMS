`locate.peaks` <-
function(peak.base, num.pts=5, R2.thresh = 0.98, oneside.min=1, 
        peak.method = "parabola", thresh=-Inf){
    loc.max <- .loc.maxes(peak.base[,2])
    loc.max <- loc.max[peak.base[loc.max,2]>=thresh]    
    if(peak.method=="parabola"){
        loc.max <- loc.max[loc.max >= oneside.min+1 & loc.max <= length(peak.base[,2])-oneside.min]
        all.peaks <- sapply(loc.max,function(x){
            locs <- x + ((-num.pts+oneside.min+1):(num.pts-oneside.min-1))
            locs <- locs[locs>=1 & locs<=dim(peak.base)[1]]
            .peak.parab(peak.base[locs,], num.pts, R2.thresh)
        })
        all.peaks <- data.frame(matrix(unlist(all.peaks),ncol=3,byrow=TRUE))
        names(all.peaks) <- c("Center_hat","Max_hat","Width_hat")     
        all.peaks <- all.peaks[all.peaks$Width_hat>0,]
        rownames(all.peaks) <- 1:dim(all.peaks)[1]
        for(i in 1:3){all.peaks[,i]<-as.numeric(all.peaks[,i])}
    } else {
        all.peaks <- data.frame(peak.base[loc.max,], NA)
        names(all.peaks) <- c("Center_hat","Max_hat","Width_hat")     
    }
    all.peaks
}

