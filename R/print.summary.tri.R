print.summary.tri<-function(summ.tri.obj)
  {
    cat("triangulation:\n")
    cat("Call:", deparse(summ.tri.obj$call),"\n")
    cat("number of nodes:",summ.tri.obj$n,"\n")
    cat("number of arcs:",summ.tri.obj$na,"\n")
    cat("number of boundary nodes:",summ.tri.obj$nb,"\n")
    cat("number of triangles:",summ.tri.obj$nt,"\n")
    cat("number of constraints:",summ.tri.obj$nc,"\n")
  }
