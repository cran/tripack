.First.lib <- function(lib, pkg) library.dynam("tripack", pkg, lib)


if(version$minor < "62") library.dynam("tripack")
