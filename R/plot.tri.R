plot.tri<-function(x,...,add=FALSE,xlim=range(x$x),
                   ylim=range(x$y),do.points=TRUE)
{
  if(!inherits(x,"tri"))
    stop("x must be of class \"tri\"")
  tnabor<- integer(x$tlnew)
  nnabs <- integer(x$n)
  nptr <- integer(x$n)
  nptr1 <- integer(x$n)
  nbnos <- integer(x$n)
  ans<-.Fortran("troutq",
                 as.integer(x$nc),
                 as.integer(x$lc),
                 as.integer(x$n),
                 as.double(x$x),
                 as.double(x$y),
                 as.integer(x$tlist),
                 as.integer(x$tlptr),
                 as.integer(x$tlend),
                 as.integer(6),
                 nnabs=as.integer(nnabs),
                 nptr=as.integer(nptr),
                 nptr1=as.integer(nptr1),
                 tnabor=as.integer(tnabor),
                 nbnos=as.integer(nbnos),
                 na=as.integer(0),
                 nb=as.integer(0),
                 nt=as.integer(0),
                 PACKAGE = "tripack")
  if(!add)
    {
      plot.new()
      plot.window(xlim=xlim,ylim=ylim,"")
    }
  for (i in 1:x$n)
    {
      inb<-ans$tnabor[ans$nptr[i]:ans$nptr1[i]]
      for (j in inb)
        lines(c(x$x[i],x$x[j]),c(x$y[i],x$y[j]), ...)
    }
  if(do.points) points(x$x,x$y)
  if(!add) title("Delaunay triangulation",deparse(substitute(x)))
}
