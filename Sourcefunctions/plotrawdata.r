plotrawdata <- function(observed)
{
    ### Note here the 'observed' is of the wide data format
    if(is.data.frame(observed)) nm <- names(observed) else{
        if(is.matrix(observed)) nm <- colnames(observed)
    }
    yl <- range(observed[,2:ncol(observed)])
    plot(observed[,1],observed[,2],xlab=nm[1],ylim=yl,ylab='concentration')
    for(i in 2:ncol(observed))
        points(observed[,1],observed[,i],col=i-1,pch=i-1)
    legend("topright", inset=c(0.05, 0.05), legend=nm,
      col=1:(i-1),pch=1:(i-1))
}
