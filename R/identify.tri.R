identify.tri<-function(tri.obj)
  {
    if(!inherits(tri.obj,"tri"))
      stop("tri.obj must be of class \"tri\"")
    labels<-paste("(",tri.obj$x,",",tri.obj$y,")", sep ="")
    identify(tri.obj$x,tri.obj$y,labels=labels)
  }
