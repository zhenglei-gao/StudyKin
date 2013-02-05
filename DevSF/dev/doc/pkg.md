Package Development
========================================================

This is an R Markdown document. Markdown is a simple formatting syntax for authoring web pages (click the **MD** toolbar button for help on Markdown).

When you click the **Knit HTML** button a web page will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:


```r
library(roxygen2)
```

```
## Loading required package: digest
```

```r
# upper directory of the package
setwd("C:/Users/z.gao/Documents/GitHub/KINGUI/KinGUII/DevSF/dev/doc/")
roxygenize("kinGUI2")
```


You can also embed plots, for example:


```r
plot(cars)
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 


