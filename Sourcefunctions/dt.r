## FUNCTIONS To plot the Curve to calculate DT90
## DFOP
fDFOP <- function(t, x) {
  ((g * exp( - k1 * t) + (1 - g) * exp( - k2 * t)) - (1 - x/100))^2
}

ffDFOP90 <- function(t){
  fDFOP(t,x=90)
}
## tseq <- seq(0,1000,by=0.1)
## plot(tseq,ff90(tseq))
## abline(h=0)