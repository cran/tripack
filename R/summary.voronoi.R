summary.voronoi<-function(voronoi.obj)
{
  if(!inherits(voronoi.obj,"voronoi"))
    stop("voronoi.obj must be of class \"voronoi\"")
  ans<-list(nn=length(voronoi.obj$x),
            nd=length(voronoi.obj$dummy.x),
            call=voronoi.obj$call)
  class(ans)<-"summary.voronoi"
  ans
}
