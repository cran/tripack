.First.lib <- function(lib, pkg) {
  if(getRversion()<="0.62")
    stop("This version for R 0.62 or later")
  library.dynam("tripack", pkg, lib)
}

