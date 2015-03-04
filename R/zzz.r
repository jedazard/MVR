##########################################################################################################################################
#     MVR
##########################################################################################################################################

.onAttach <- function(libname, pkgname) {
    SSver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), 
                      fields="Version")
    packageStartupMessage(paste(pkgname, SSver))
    packageStartupMessage("Type MVR.news() to see new features, changes, and bug fixes\n")
}
