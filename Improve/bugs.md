Bugs and Fixes
========================================================
In case the library cannot be loaded.

```r
library("KineticEval", lib.loc = "M:/Personal Data/R/win-library/2.15/")
library("deSolve", lib.loc = "M:/Personal Data/R/win-library/2.15/")
library("rootSolve", lib.loc = "M:/Personal Data/R/win-library/2.15/")
library("minpack.lm", lib.loc = "M:/Personal Data/R/win-library/2.15/")
library("coda", lib.loc = "M:/Personal Data/R/win-library/2.15/")
library("FME", lib.loc = "M:/Personal Data/R/win-library/2.15/")
library("mkin", lib.loc = "M:/Personal Data/R/win-library/2.15/")
library("numDeriv", lib.loc = "M:/Personal Data/R/win-library/2.15/")
library("optimx", lib.loc = "M:/Personal Data/R/win-library/2.15/")
Sys.getenv()
```


Now assuming the sourcefunctions/packages are loaded.


## Bug 1: Cannot reach the optimal solution


```r
##
mkinmodini <- mkinmod.full(Parent = list(type = "SFO", to = "Metab"), Metab = list(type = "SFO", 
    M0 = list(ini = 0, fixed = 0, lower = 0, upper = Inf)), data = andrew)
```

```
## Error: could not find function "mkinmod.full"
```

```r
Fit <- mkinfit.full(mkinmodini, ctr = kingui.control(method = "L-BFGS-B", submethod = "Port"))
```

```
## Error: could not find function "mkinfit.full"
```

```r
compare_multi_kinmod(mkinmodini)
```

```
## Error: could not find function "compare_multi_kinmod"
```

```r

compare_multi_kinmod(mkinmodini, rbind(Fit$par, c(2.0633, 0.2105, 0.3033, 0.8945, 
    1)))
```

```
## Error: could not find function "compare_multi_kinmod"
```

```r
newmod1 <- update_kinmod(mkinmodini, newparms = Fit$par)$new_mkinmod
```

```
## Error: could not find function "update_kinmod"
```

```r
Fit1 <- mkinfit.full(newmod1, ctr = kingui.control(method = "L-BFGS-B", submethod = "Port"))
```

```
## Error: could not find function "mkinfit.full"
```

```r
newmod2 <- update_kinmod(mkinmodini, newparms = c(2.0633, 0.2105, 0.3033, 0.8945, 
    1))$new_mkinmod
```

```
## Error: could not find function "update_kinmod"
```

```r
Fit2 <- mkinfit.full(newmod2, ctr = kingui.control(method = "L-BFGS-B", submethod = "Port"))
```

```
## Error: could not find function "mkinfit.full"
```

```r
compare_multi_kinmod(mkinmodini, rbind(Fit$par, c(2.0633, 0.2105, 0.3033, 0.8945, 
    1), Fit2$par))
```

```
## Error: could not find function "compare_multi_kinmod"
```


## Bug 2: the optimized parameters are not the same as specified in the model

```r
## a possible bug that the optimized parameters are not the same as
## specified in the model
mkinmodini <- mkinmod.full(Parent = list(type = "SFO", to = "Metab"), Metab = list(type = "SFO"), 
    data = andrew)
```

```
## Error: could not find function "mkinmod.full"
```

```r
Fit <- mkinfit.full(mkinmodini, ctr = kingui.control(method = "L-BFGS-B", submethod = "Port"))
```

```
## Error: could not find function "mkinfit.full"
```

```r
newmod1 <- update_kinmod(mkinmodini, newparms = Fit$par)$new_mkinmod
```

```
## Error: could not find function "update_kinmod"
```

```r
Fit1 <- mkinfit.full(newmod1, ctr = kingui.control(method = "L-BFGS-B", submethod = "Port"))
```

```
## Error: could not find function "mkinfit.full"
```

```r
newmod2 <- update_kinmod(mkinmodini, newparms = c(2.0633, 0.2105, 0.3033, 0.8945, 
    1))$new_mkinmod
```

```
## Error: could not find function "update_kinmod"
```

```r
Fit2 <- mkinfit.full(newmod2, ctr = kingui.control(method = "L-BFGS-B", submethod = "Port"))
```

```
## Error: could not find function "mkinfit.full"
```

```r
compare_multi_kinmod(mkinmodini, rbind(Fit$par, c(2.0633, 0.2105, 0.3033, 0.8945, 
    1), Fit2$par))
```

```
## Error: could not find function "compare_multi_kinmod"
```

You can also embed plots, for example:




