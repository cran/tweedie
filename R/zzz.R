# R/zzz.R

.onLoad <- function(libname, pkgname) {
  # This explicitly ensures the shared object is loaded for the current namespace
  # This is a safer way to get the DLLInfo object R is complaining about.
#  utils::getFromNamespace("loadNamespace", "base")(pkgname)
}

.onUnload <- function(libpath) {
  # library.dynam.unload uses the package name
  library.dynam.unload("tweedie", libpath)
}