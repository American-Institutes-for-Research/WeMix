# @author Paul Bailey
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("WeMix v", utils::packageDescription("WeMix")$Version, "\n"))
}
