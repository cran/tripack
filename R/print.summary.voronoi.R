print.summary.voronoi<-function(summ.vo.obj)
  {
    cat("voronoi mosaic\n")
    cat("Call:", deparse(summ.vo.obj$call),"\n")
    cat(summ.vo.obj$nn, "nodes\n")
    cat(summ.vo.obj$nd, "dummy nodes\n")
  }
