plot.tri<-function(tri.obj,add=F,xlim=range(tri.obj$x),
                   ylim=range(tri.obj$y),do.points=T,...)
{
  if(!inherits(tri.obj,"tri"))
    stop("tri.obj must be of class \"tri\"")
  tnabor<-matrix(-1,tri.obj$n,tri.obj$n)
  storage.mode(tnabor)<-"integer"
  ans<-.Fortran("troutp",
                 as.integer(tri.obj$nc),
                 as.integer(tri.obj$lc),
                 as.integer(tri.obj$n),
                 as.double(tri.obj$x),
                 as.double(tri.obj$y),
                 as.integer(tri.obj$tlist),
                 as.integer(tri.obj$tlptr),
                 as.integer(tri.obj$tlend),
                 as.integer(6),
                 tnabor=as.integer(tnabor),
                 na=as.integer(0),
                 nb=as.integer(0),
                 nt=as.integer(0)
                 )
  ans$tnabor<-matrix(ans$tnabor,tri.obj$n,tri.obj$n)
  if(!add)
    {
      plot.new()
      plot.window(xlim=xlim,ylim=ylim,"")
    }
  for (i in 1:tri.obj$n)
    {
      inb<-ans$tnabor[i,ans$tnabor[i,]>0]
      for (j in inb)
        lines(c(tri.obj$x[i],tri.obj$x[j]),c(tri.obj$y[i],tri.obj$y[j]), ...)
    }
  if(do.points) points(tri.obj$x,tri.obj$y)
  if(!add) title("Delaunay triangulation",deparse(substitute(tri.obj)))
}
