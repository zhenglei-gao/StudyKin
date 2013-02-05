setwd('C:/Users/gbbfx/Work/Projects/Kinetics/KINGUI/package/workingcopy/DevSF/dev/doc')
require( highlight )
driver <- HighlightWeaveLatex(box=TRUE)
Sweave( 'kinGUI2.Rnw', driver = driver )
###################################
library(Rd2roxygen)
options(roxygen.comment = "##' ")
(info = parse_file('C:/Users/gbbfx/Work/Projects/Kinetics/KINGUI/package/workingcopy/DevSF/dev/man/mkinfit.gui.Rd'))
cat(create_roxygen(info), sep = "\n") # parse_and_save() combines these two steps
(info = parse_file('C:/Users/gbbfx/Work/Projects/Kinetics/KINGUI/package/workingcopy/DevSF/dev/man/summary.kingui0.Rd'))
cat(create_roxygen(info), sep = "\n")
(info = parse_file('C:/Users/gbbfx/Work/Projects/Kinetics/KINGUI/package/workingcopy/DevSF/dev/man/IRLSkinfit.gui.Rd'))
cat(create_roxygen(info), sep = "\n")
(info = parse_file('C:/Users/gbbfx/Work/Projects/Kinetics/KINGUI/package/workingcopy/DevSF/dev/man/mcmckinfit.gui.Rd'))
cat(create_roxygen(info), sep = "\n")
(info = parse_file('C:/Users/gbbfx/Work/Projects/Kinetics/KINGUI/package/workingcopy/DevSF/dev/man/plot.mcmckingui0.Rd'))
cat(create_roxygen(info), sep = "\n")
library(roxygen2)
(info = parse_file('C:/Users/gbbfx/Work/Projects/Kinetics/KINGUI/package/workingcopy/DevSF/dev/man/kingui.control.Rd'))
cat(create_roxygen(info), sep = "\n")
library(roxygen2)
(info = parse_file('C:/Users/gbbfx/Work/Projects/Kinetics/KINGUI/package/workingcopy/DevSF/dev/man/mkinmod.gui.Rd'))
cat(create_roxygen(info), sep = "\n")
library(roxygen2)
(info = parse_file('C:/Users/gbbfx/Work/Projects/Kinetics/KINGUI/package/workingcopy/DevSF/dev/man/kingraph.Rd'))
cat(create_roxygen(info), sep = "\n")
library(roxygen2)
## example usage ##package.skeleton("pkg")
package.skeleton('kinGUI2',code_files=c('C:/Users/gbbfx/Work/Projects/Kinetics/KINGUI/package/workingcopy/Sourcefunctions/IRLSkinfit.full.R','C:/Users/gbbfx/Work/Projects/Kinetics/KINGUI/package/workingcopy/Sourcefunctions/IRLSkinfit.gui.r','C:/Users/gbbfx/Work/Projects/Kinetics/KINGUI/package/workingcopy/Sourcefunctions/KinGUI2-package.r','C:/Users/gbbfx/Work/Projects/Kinetics/KINGUI/package/workingcopy/Sourcefunctions/mkinfit.gui.R','C:/Users/gbbfx/Work/Projects/Kinetics/KINGUI/package/workingcopy/Sourcefunctions/mkinfit.full.R','C:/Users/gbbfx/Work/Projects/Kinetics/KINGUI/package/workingcopy/Sourcefunctions/mcmckinfit.full.r','C:/Users/gbbfx/Work/Projects/Kinetics/KINGUI/package/workingcopy/Sourcefunctions/kinplot.R','C:/Users/gbbfx/Work/Projects/Kinetics/KINGUI/package/workingcopy/Sourcefunctions/kingraph.R','C:/Users/gbbfx/Work/Projects/Kinetics/KINGUI/package/workingcopy/Sourcefunctions/mkinmod.full.r','C:/Users/gbbfx/Work/Projects/Kinetics/KINGUI/package/workingcopy/Sourcefunctions/kingui.control.r'),force=TRUE)
roxygenize('kinGUI2')
## R CMD Rd2dvi --pdf --title="Package KinGUI2" -o KinGUI2_help.pdf man/*.Rd
###################################
## Read './pkg/Read-and-delete-me' file, compile the DESCRIPTION fiels according to your needs and delete './pkg/Read-and-delete-me'.
## Now the devtools magic:
library("devtools")
pkg <- as.package("KinGUI2")

load_all(pkg, reset=T) # to reload the package without having to restart R
document(pkg) # to be used together with roxygen2 to creating the corresponding Rd files
run_examples(pkg) # to check the examples for the different functions
devtools:::check(pkg) # to verified if your package raises errors or warnings
devtools:::build(pkg)

## install(pkg) # install your package
# release()

## from formatR package
tidy.source = function(file = choose.files()) {
   exprs = parse(file)
   for (i in 1:length(exprs)) {
       dep = paste(deparse(exprs[i]), collapse = "\n")
       dep = substring(dep, 12, nchar(dep) - 1)
       cat(dep, "\n")
   }
}
### Example of generating calling graphs #####
library(mvbutils)
inner1 <- function() {
  "This is the inner1 function"
}

inner2 <- function() {
  "This is the inner2 function"
}

outer <- function(x) {
   i1 <- inner1()
   i2 <- inner2()
}
foodweb()
help(foodweb)
foodweb(border = TRUE, expand.xbox = 3,
        boxcolor = "#FC6512", textcolor = "black",
        cex = 1.2, lwd=2)
foodweb(where = "package:survival", prune = "survexp",
        border = TRUE,
        expand.xbox = 2, boxcolor = "#FC6512",
        textcolor = "black", cex = 1.0, lwd=2)
mtext("The survexp function foodweb")
########################################################

sourceDir('S:/Software & EDV/Models/FOCUS Kinetics/KinGUII/R-2.12.1/SourceFunctions/')
foodweb(c('IRLSkinfit.full','mkinfit.full','mcmckinfit.full','modFit1','modMCMC','mkinmod.full','completeCompound', "kingui.control"  , "ForwardCalcFF" , "modCost" ,'summary.kingui','summary.mcmckingui','kingraph','chi2err'),ancestors=FALSE,descendents=FALSE, border =FALSE,expand.xbox = 1.2, boxcolor = "grey",textcolor = "black", cex = 0.9, lwd=2)
mtext(expression(paste("Calling Graph for package ", bold(KinGUII),sep=' ')))
dev.print(pdf,'callgraph.KinGUII.pdf')
foodweb( c( find.funs(),find.funs("package:FME"),find.funs('package:mkin')), prune = "mkinmod.full",
        border = TRUE,
        expand.xbox = 2, boxcolor = "#FC6512",
        textcolor = "black", cex = 1.0, lwd=2)
mtext(expression(paste("Calling Graph for function ", bold(mkinmod.full),sep=' ')))
dev.print(pdf,'callgraph.mkinmod.full.pdf')
###########################################################
foodweb(c( find.funs(),find.funs("package:FME"),'solnp'),where=1, prune = "IRLSkinfit.full",border = TRUE,expand.xbox = 2, boxcolor = "#FC6512",textcolor = "black", cex = 1.0, lwd=2)
mtext(expression(paste("Calling Graph for function ", bold(IRLSkinfit.full),sep=' ')))
dev.print(pdf,'call.IRLSkinfit.pdf')

foodweb(c( find.funs(),find.funs("package:FME"),'solnp'), prune = "mcmckinfit.full",
        border = TRUE,
        expand.xbox = 1.2, boxcolor = "grey",
        textcolor = "black", cex = 0.9, lwd=2)
mtext(expression(paste("Calling Graph for function ", bold(mcmckinfit.full),sep=' ')))
dev.print(pdf,'call.mcmckinfit.pdf')

foodweb(c( find.funs(),find.funs("package:FME"),'solnp','optimx','nls.lm','nlm','nlminb','bobyqa'), prune = "modFit1",ancestors=FALSE,
        border = FALSE,
        expand.xbox = 1.5, boxcolor = "grey",
        textcolor = "black", cex = 0.9, lwd=2)
mtext(expression(paste("Calling Graph for function ", bold(modFit1),sep=' ')))
dev.print(pdf,'call.modFit1.pdf')

foodweb(c( find.funs(),'solnp'), prune = "mkinfit.full",
        border = TRUE,
        expand.xbox = 2, boxcolor = "#FC6512",
        textcolor = "black", cex = 1.0, lwd=2)
mtext(expression(paste("Calling Graph for function ", bold(mkinfit.full),sep=' ')))
dev.print(pdf,'call.mkinfit.pdf')



foodweb(c('summary.kingui','print.summary.kingui','myformat'), prune = "summary.kingui",
        border = TRUE,
        expand.xbox = 2, boxcolor = "#FC6512",
        textcolor = "black", cex = 1.0, lwd=2)
mtext(expression(paste("Calling Graph for function ", bold(summary.kingui),sep=' ')))
dev.print(pdf,'call.summary.kingui.pdf')
##
foodweb(c('summary.mcmckingui','print.summary.mcmckingui','myformat' ,find.funs("package:FME"),'mcmc'), prune = "summary.mcmckingui",
        border = TRUE,
        expand.xbox = 1.5, boxcolor = "#FC6512",
        textcolor = "black", cex = 1.0, lwd=2)
mtext(expression(paste("Calling Graph for function ", bold(summary.mcmckingui),sep=' ')))
dev.print(pdf,'call.summary.mcmckingui.pdf')

##################################################
foodweb(c( find.funs(),find.funs("package:FME"),'solnp'), prune = "modCost",
        border = TRUE,
        expand.xbox = 2, boxcolor = "#FC6512",
        textcolor = "black", cex = 1.0, lwd=2)
mtext(expression(paste("Calling Graph for function ", bold(modCost),sep=' ')))
#########
require(tikzDevice)
require(Hmisc)
