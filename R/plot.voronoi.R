"plot.voronoi" <- function(v.obj,add=F,
                           xlim=c(min(v.obj$tri$x)-
                             0.1*diff(range(v.obj$tri$x)),
                             max(v.obj$tri$x)+
                             0.1*diff(range(v.obj$tri$x))),
                           ylim=c(min(v.obj$tri$y)-
                             0.1*diff(range(v.obj$tri$y)),
                             max(v.obj$tri$y)+
                             0.1*diff(range(v.obj$tri$y))),
                           all=F,
                           do.points=T,
                           main="Voronoi mosaic",
                           sub=deparse(substitute(v.obj)),
                           ...)
  {

    
    if(all)
      {
        xlim<-c(min(v.obj$x)-0.1*diff(range(v.obj$x)),
                max(v.obj$x)+0.1*diff(range(v.obj$x)))
        ylim<-c(min(v.obj$y)-0.1*diff(range(v.obj$y)),
                max(v.obj$y)+0.1*diff(range(v.obj$y)))
      }
    
    n<-length(v.obj$x)

    if(!add)
      {
        plot.new()
        plot.window(xlim=xlim,ylim=ylim,"")
      }

    if(do.points) points(v.obj$x,v.obj$y)

    for (i in 1:n)
      {
        if(v.obj$node[i])
          # Triangle i has positive area.
          # Connect circumcircle center of triangle i with neighbours:
          {
            # Find neighbour triangles
            tns<-sort(c(v.obj$n1[i],v.obj$n2[i],v.obj$n3[i]))
            for(j in 1:3)
              {
                # Connect (if triangle exists and has positive area).
                if(tns[j]>0)
                  {
                  # simple node
                    if(v.obj$node[tns[j]])
                      lines(c(v.obj$x[i],v.obj$x[tns[j]]),
                            c(v.obj$y[i],v.obj$y[tns[j]]),...)
                  }
                else if(tns[j]<0){
                  # dummy node
                  lines(c(v.obj$x[i],v.obj$dummy.x[-tns[j]]),
                        c(v.obj$y[i],v.obj$dummy.y[-tns[j]]),
                        lty="dashed",...) }
              }
          }
      }
    if(!add)
      title(main = main, sub =sub)
  }
