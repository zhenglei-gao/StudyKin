if ("devtools" %in% loadedNamespaces()) {
  
  stop("You must restart R before installing devtools")
  
}



url <- "https://gist.github.com/raw/4506250/devtools.zip"

temp <- file.path(tempdir(), "devtools.zip")



setInternet2(TRUE)

suppressWarnings(download.file(url, temp, mode = "wb"))

install.packages(temp, repos = NULL)

file.remove(temp)

########################
library(devtools)

#########################
library(roxygen2)
#########################
Install the package