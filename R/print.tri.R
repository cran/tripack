print.tri<-function(tri.obj)
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
  cat("triangulation nodes with neigbours:\n")
  cat("node: (x,y): neighbours\n")
  for (i in 1:tri.obj$n)
    {
      cat(i,": (",tri.obj$x[i],",",tri.obj$y[i],") [",ans$nnabs[i],"]",sep="")
      cat(":",sort(ans$tnabor[ans$nptr[i]:ans$nptr1[i]]),"\n",sep=" ")
    }
  cat("number of nodes:",tri.obj$n,"\n")
  cat("number of arcs:",ans$na,"\n")
  cat("number of boundary nodes:",ans$nb,"\n")
  cat("boundary nodes: ",ans$nbnos[1:ans$nb], "\n", sep=" ")
  cat("number of triangles:",ans$nt,"\n")
  cat("number of constraints:",tri.obj$nc,"\n")
}
