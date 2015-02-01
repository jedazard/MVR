##########################################################################################################################################
#     MVR
##########################################################################################################################################

.onAttach <- function(libname, pkgname) {
    SSver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), 
                      fields="Version")
    packageStartupMessage("\n")
    packageStartupMessage(paste(pkgname, SSver))
    packageStartupMessage("\n")
    packageStartupMessage("\n")
    packageStartupMessage("Type MVR.news() to see new features, changes, and bug fixes.")
    packageStartupMessage("\n")
    packageStartupMessage("\n")
}
