.onAttach <- function(...) {

    packageStartupMessage("\n## ----------------------------------------- ##")
    packageStartupMessage("## L1 Norm Ideal Point Estimation Package")
    packageStartupMessage("## Sooahn Shin, Johan Lim, and Jong Hee Park")
    packageStartupMessage("## ----------------------------------------- ##")
    
}

.onUnload <- function(libpath) {
    library.dynam.unload("l1ideal", libpath)
}
