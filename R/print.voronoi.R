print.voronoi<-function(voronoi.obj)
{
  if(!inherits(voronoi.obj,"voronoi"))
    stop("voronoi.obj must be of class \"voronoi\"")
  cat("voronoi mosaic:\n")
  cat("nodes: (x,y): neighbours (<0: dummy node)\n")
  for (i in 1:length(voronoi.obj$x))
    {
      if(voronoi.obj$node[i]){
        cat(i,": (",voronoi.obj$x[i],",",voronoi.obj$y[i],")",sep="")
        cat(":",voronoi.obj$n1[i],voronoi.obj$n2[i],voronoi.obj$n3[i],"\n",sep=" ")
      }
    }
  cat("dummy nodes: (x,y)\n")
  for (i in 1:length(voronoi.obj$dummy.x))
    {
      cat(i,": (",voronoi.obj$dummy.x[i],",",voronoi.obj$dummy.y[i],")\n",sep="")
    }

}
