# CHANGES IN KineticEval VERSION 0.10

## NEW FEATURES

- a new function `KineticEval()` which is similar to **brew** and mustache, e.g. it expands `pi is {{pi}}` to `pi is 3.14`; it can also be used for building child documents (see https://github.com/yihui/knitr-examples/blob/master/075-knit-expand.Rnw for example) (#397) (thanks, Frank Harrell)

- `mkinmod.full()` gained a new argument `encoding` to blabla

- a new function `Sweave2knitr()` to convert Sweave documents to **knitr**; several automatic translations can be done, e.g. `results=tex` to `results='asis'`, `width=5` to `fig.width=5`, `echo=true` to `echo=TRUE`, `keep.source=TRUE` to `tidy=FALSE`, `eps=TRUE` to `dev='postscript'`, `\SweaveOpts{...}` to `opts_chunk$set(...)` and so on; see the documentation in the package for details (#451)

- the chunk label is used as the id of the div element in R HTML output, e.g. `<div id='chunk-label'>...</div>`

## MAJOR CHANGES

- (IMPORTANT) Metablate with a **DFOP** degradation model was ...

- accordingly, the pattern elements `global.options` and `inline.doc` were removed from `knit_patterns` (`\SweaveOpts{}` and `\SweaveInput{}` will no longer be supported; please call `Sweave2knitr()` to convert incompatible Sweave documents)


## MINOR CHANGES

- for inline R code, the value is returned only if the R code prints a visible value, e.g. `\Sexpr{x <- 1}` will be empty, and `\Sexpr{pi}` will return the value of pi

## BUG FIXES

- fixed #1: no longer uses `\\\\` in LaTeX output; only a single line break is converted to `\\` (thanks, Kevin Wright)

- `render_html()` guarantees that the R source code is highlighted when the chunk option `highlight = TRUE` (#447) (thanks, Ramnath Vaidyanathan)

- `dep_auto()` was unable to find the cache files if the input document is not under the current working directory (thanks, Hui Yao)


  
# CHANGES IN KineticEval VERSION 0.1

## NEW FEATURES
  	
- first version of KineticEval: it covers most features in **KinGUI** and **mkinfit**; see package homepage for documentation and examples: http://KinGUII/

## MISC

- **KinGUII** was distributed via DropBox before.

- in this NEWS file, #n means the issue number on GitHub, e.g. #11 is https://github.com/zhenglei-gao/KinGUII/issues/11

