plot.tri<-function(tri.obj,add=F,xlim=range(tri.obj$x),
                   ylim=range(tri.obj$y),do.points=T,...)
{
  if(!inherits(tri.obj,"tri"))
    stop("tri.obj must be of class \"tri\"")
  tnabor<- integer(tri.obj$tlnew)
  nnabs <- integer(tri.obj$n)
  nptr <- integer(tri.obj$n)
  nptr1 <- integer(tri.obj$n)
  nbnos <- integer(tri.obj$n)
  ans<-.Fortran("troutq",
                 as.integer(tri.obj$nc),
                 as.integer(tri.obj$lc),
                 as.integer(tri.obj$n),
                 as.double(tri.obj$x),
                 as.double(tri.obj$y),
                 as.integer(tri.obj$tlist),
                 as.integer(tri.obj$tlptr),
                 as.integer(tri.obj$tlend),
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
  for (i in 1:tri.obj$n)
    {
      inb<-ans$tnabor[ans$nptr[i]:ans$nptr1[i]]
      for (j in inb)
        lines(c(tri.obj$x[i],tri.obj$x[j]),c(tri.obj$y[i],tri.obj$y[j]), ...)
    }
  if(do.points) points(tri.obj$x,tri.obj$y)
  if(!add) title("Delaunay triangulation",deparse(substitute(tri.obj)))
}
