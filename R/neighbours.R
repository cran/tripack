neighbours<-function(tri.obj)
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
  ret<-rep(NULL,tri.obj$n)
  for (i in 1:tri.obj$n)
    {
#      cat(i,": (",tri.obj$x[i],",",tri.obj$y[i],")",sep="")
#      cat(":",ans$tnabor[i,ans$tnabor[i,]!=-1],"\n",sep=" ")
      ret[i]<-ans$tnabor[i,ans$tnabor[i,]!=-1]
    }
  ret
}
