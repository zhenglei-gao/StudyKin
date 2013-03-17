Trust Region Reflective?
========================================================




Other's implementation

* https://bitbucket.org/carandraug/octave/src/13d1e9bfa362/scripts/optimization/lsqnonneg.m?at=default

-------------------------

Stackoverflow Question
===========================================================
My original question was "Implementation of trust-region-reflective optimization algorithm in R". However, on the way of producing a reproducible example(thanks @Ben for his advice), I realize that my problem is that in Matlab, one function `lsqnonlin` is good(meaning no need to choose a good starting value, fast) enough for most cases I have, while in R, there is not such a one-for-all function. Different optmization algorithm works well in different cases. Different algorithms reach different solutions. The reason behind this may not be that the optimization algorithms in R is inferior to the trust-region-reflective algorithm in Matlab, it could also be related to how R handles Automatic Differentiation. This problem comes actually from interrupted work two years ago. At that time, Prof. John C Nash, one of the authors of the package **optimx** has suggested that there has been quite a lot of work for Matlab for Automatic Differentiation, which might be the reason that the Matlab lsqnonlin performs better than the optimization functions/algorithms in R. I am not able to figure it out with my knowledge. 

The example below shows some problems I have encountered(More reproducible examples are coming). To run the examples, first run `install_github("KineticEval","zhenglei-gao")`. You need to install package **mkin** and its dependencies and may also need to install a bunch of other packages for different optimization algorithms.

Basically I am trying to solve nonlinear least-squares curve fitting problems as described in the Matlab function `lsqnonlin` 's documentation (http://www.mathworks.de/de/help/optim/ug/lsqnonlin.html). The curves in my case are modeled by a set of differential equations. I will explain a bit more with the examples. Optimization algorithms I have tried including:

* Marq from `nls.lm`, the Levenburg-Marquardt
* Port from `nlm.inb`
* L-BGFS-B from `optim`
* spg from `optimx`
* `solnp` of package **Rsolnp**

I have also tried a few others but here are the ones I will often use. 

## Example 1: A simple case

I will give the R codes first and explain later.


```r
ex1 <- mkinmod.full(Parent = list(type = "SFO", to = "Metab", sink = TRUE, k = list(ini = 0.1, 
    fixed = 0, lower = 0, upper = Inf), M0 = list(ini = 195, fixed = 0, lower = 0, 
    upper = Inf), FF = list(ini = c(0.1), fixed = c(0), lower = c(0), upper = c(1)), 
    time = c(0, 2.8, 6.2, 12, 29.2, 66.8, 99.8, 127.5, 154.4, 229.9, 272.3, 
        288.1, 322.9), residue = c(157.3, 206.3, 181.4, 223, 163.2, 144.7, 85, 
        76.5, 76.4, 51.5, 45.5, 47.3, 42.7), weight = c(1, 1, 1, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 1)), Metab = list(type = "SFO", k = list(ini = 0.1, fixed = 0, 
    lower = 0, upper = Inf), M0 = list(ini = 0, fixed = 1, lower = 0, upper = Inf), 
    residue = c(0, 0, 0, 1.6, 4, 12.3, 13.5, 12.7, 11.4, 11.6, 10.9, 9.5, 7.6), 
    weight = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)))
ex1$diffs
Fit <- NULL
alglist <- c("L-BFGS-B", "Marq", "Port", "spg", "solnp")
for (i in 1:5) {
    Fit[[i]] <- mkinfit.full(ex1, plot = TRUE, quiet = TRUE, ctr = kingui.control(method = alglist[i], 
        submethod = "Port", maxIter = 100, tolerance = 1e-06, odesolver = "lsoda"))
}
names(Fit) <- alglist
(lapply(Fit, function(x) x$par))
kinplot(Fit[[2]])
unlist(lapply(Fit, function(x) x$ssr))
```

The output from the last line is:

    L-BFGS-B     Marq     Port      spg    solnp 
    5735.744 4714.500 5780.446 5728.361 4714.499 

Except for "Marq" and "solnp", the other algorithms did not reach the optimum.Besides, 'spg' method (also other methods like 'bobyqa') need too many function evaluations for such a simple case . Moreover, if I change the starting value and make `k_Parent=0.0058` (the optimum value for that parameter) instead of the random choosen `0.1`,  "Marq" cannot find the optimum any more! (Code provided below). I have also had datasets where "solnp" does not find the optimum. However, if I use `lsqnonlin` in Matlab, I haven't encountered any difficulties for such simple cases. 


```r
ex1_a <- mkinmod.full(Parent = list(type = "SFO", to = "Metab", sink = TRUE, 
    k = list(ini = 0.0058, fixed = 0, lower = 0, upper = Inf), M0 = list(ini = 195, 
        fixed = 0, lower = 0, upper = Inf), FF = list(ini = c(0.1), fixed = c(0), 
        lower = c(0), upper = c(1)), time = c(0, 2.8, 6.2, 12, 29.2, 66.8, 99.8, 
        127.5, 154.4, 229.9, 272.3, 288.1, 322.9), residue = c(157.3, 206.3, 
        181.4, 223, 163.2, 144.7, 85, 76.5, 76.4, 51.5, 45.5, 47.3, 42.7), weight = c(1, 
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)), Metab = list(type = "SFO", k = list(ini = 0.1, 
    fixed = 0, lower = 0, upper = Inf), M0 = list(ini = 0, fixed = 1, lower = 0, 
    upper = Inf), residue = c(0, 0, 0, 1.6, 4, 12.3, 13.5, 12.7, 11.4, 11.6, 
    10.9, 9.5, 7.6), weight = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)))

Fit_a <- NULL
alglist <- c("L-BFGS-B", "Marq", "Port", "spg", "solnp")
for (i in 1:5) {
    Fit_a[[i]] <- mkinfit.full(ex1_a, plot = TRUE, quiet = TRUE, ctr = kingui.control(method = alglist[i], 
        submethod = "Port", maxIter = 100, tolerance = 1e-06, odesolver = "lsoda"))
}
names(Fit_a) <- alglist
lapply(Fit_a, function(x) x$par)
unlist(lapply(Fit_a, function(x) x$ssr))
```

Now the output from last line is:

    L-BFGS-B     Marq     Port      spg    solnp 
    5653.132 4866.961 5653.070 5635.372 4714.499 

I will explain what I am optimising here. If you have run the above script and see the curves, we use a two-compartment model with first order reactions to describe the curves. The differential equations to express the model are:


```r
ex1$diffs
```

                                                                 Parent 
                                        "d_Parent = - k_Parent * Parent" 
                                                                   Metab 
    "d_Metab = - k_Metab * Metab + k_Parent * f_Parent_to_Metab * Parent" 

For this simple case, from the differential equations we can derive the equations to describe the two curves. The to be optimized parameters are $M_0,k_p, k_m, c=\mbox{FF_parent_to_Met} $ with the constraints $M_0>0,k_p>0, k_m>0, 1> c >0$.

$$
\begin{split}
            y_{1j}&= M_0e^{-k_pt_i}+\epsilon_{1j}\\
            y_{2j} &= cM_0k_p\frac{e^{-k_mt_i}-e^{-k_pt_i}}{k_p-k_m}+\epsilon_{2j}
            \end{split}
$$

Therefore we can fit the curve without solving differential equations. 


```r
BCS1.l <- mkin_wide_to_long(BCS1)
BCS1.l <- na.omit(BCS1.l)
indi <- c(rep(1, sum(BCS1.l$name == "Parent")), rep(0, sum(BCS1.l$name == "Metab")))
sysequ.indi <- function(t, indi, M0, kp, km, C) {
    y <- indi * M0 * exp(-kp * t) + (1 - indi) * C * M0 * kp/(kp - km) * (exp(-km * 
        t) - exp(-kp * t))
    y
}
M00 <- 100
kp0 <- 0.1
km0 <- 0.01
C0 <- 0.1
library(nlme)
result1 <- gnls(value ~ sysequ.indi(time, indi, M0, kp, km, C), data = BCS1.l, 
    start = list(M0 = M00, kp = kp0, km = km0, C = C0), control = gnlsControl())
# result3 <- gnls(value ~
# sysequ.indi(time,indi,M0,kp,km,C),data=BCS1.l,start=list(M0=M00,kp=kp0,km=km0,C=C0),weights
# = varIdent(form=~1|name)) Coefficients: M0 kp km C 1.946170e+02
# 5.800074e-03 8.404269e-03 2.208788e-01
```


Doing this way, the elapsed time is almost 0, and the optimum is reached. However, we do not always have this simple case. The model can be complex and solving the differential equations is needed. See example 2

## Example 2, a complex model

I worked on this dataset a long time ago and haven't got time to finish running the following script myself. (You might need hours to finish running.) 


```r
data(BCS2)
ex2 <- mkinmod.full(Parent = list(type = "SFO", to = c("Met1", "Met2", "Met4", 
    "Met5"), k = list(ini = 0.1, fixed = 0, lower = 0, upper = Inf), M0 = list(ini = 100, 
    fixed = 0, lower = 0, upper = Inf), FF = list(ini = c(0.1, 0.1, 0.1, 0.1), 
    fixed = c(0, 0, 0, 0), lower = c(0, 0, 0, 0), upper = c(1, 1, 1, 1))), Met1 = list(type = "SFO", 
    to = c("Met3", "Met4")), Met2 = list(type = "SFO", to = c("Met3")), Met3 = list(type = "SFO"), 
    Met4 = list(type = "SFO", to = c("Met5")), Met5 = list(type = "SFO"), data = BCS2)
ex2$diffs
Fit2 <- NULL
alglist <- c("L-BFGS-B", "Marq", "Port", "spg", "solnp")
for (i in 1:5) {
    Fit2[[i]] <- mkinfit.full(ex2, plot = TRUE, quiet = TRUE, ctr = kingui.control(method = alglist[i], 
        submethod = "Port", maxIter = 100, tolerance = 1e-06, odesolver = "lsoda"))
}
kinplot(Fit2[[5]])
names(Fit) <- alglist
(lapply(Fit, function(x) x$par))
unlist(lapply(Fit, function(x) x$ssr))

```


This is an example where you will see warning messages like:

    DLSODA-  At T (=R1) and step size H (=R2), the    
      corrector convergence failed repeatedly     
      or with ABS(H) = HMIN   
    In above message, R = 
    [1] 0.000000e+00 2.289412e-09

## Issues in R

* Different algorithms reach different solutions. 
* Missing global minimum
* Lower Efficiency

I have to adimit that this problem comes actually from interrupted work two years ago. At that time, Prof. John C Nash, one of the authors of the package `optimx` has suggested that there has been quite a lot of work for Matlab for Automatic Differentiation, which might be the reason that the Matlab `lsqnonlin` performs better than the optimization functions/algorithms in R. 

Thanks @Ben Bolke at Stackoverflow for his advice of giving a **concrete, reproducible** example of a problem where R's optimizers do significantly worse than `lsqnonlin` and reframe my original question of "can I implement algorithm xxyy in R?". This makes me sit down and digging out emails and examples and have a review of my previous work. 


```r
summary(cars)
```

```
##      speed           dist    
##  Min.   : 4.0   Min.   :  2  
##  1st Qu.:12.0   1st Qu.: 26  
##  Median :15.0   Median : 36  
##  Mean   :15.4   Mean   : 43  
##  3rd Qu.:19.0   3rd Qu.: 56  
##  Max.   :25.0   Max.   :120
```


You can also embed plots, for example:


```r
###
```


## Sources

* http://stackoverflow.com/questions/13069627/library-for-trust-region-reflective-algorithms-in-c
* Ipopt
* http://www.mathworks.de/de/help/optim/ug/constrained-nonlinear-optimization-algorithms.html
* http://stackoverflow.com/questions/5527145/convert-matlab-code-to-r
