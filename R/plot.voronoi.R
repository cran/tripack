"plot.voronoi" <- function(x,...,add=FALSE,
                           xlim=c(min(x$tri$x)-
                             0.1*diff(range(x$tri$x)),
                             max(x$tri$x)+
                             0.1*diff(range(x$tri$x))),
                           ylim=c(min(x$tri$y)-
                             0.1*diff(range(x$tri$y)),
                             max(x$tri$y)+
                             0.1*diff(range(x$tri$y))),
                           all=FALSE,
                           do.points=TRUE,
                           main="Voronoi mosaic",
                           sub=deparse(substitute(x))
                           )
  {

    
    if(all)
      {
        xlim<-c(min(x$x)-0.1*diff(range(x$x)),
                max(x$x)+0.1*diff(range(x$x)))
        ylim<-c(min(x$y)-0.1*diff(range(x$y)),
                max(x$y)+0.1*diff(range(x$y)))
      }
    
    n<-length(x$x)

    if(!add)
      {
        plot.new()
        plot.window(xlim=xlim,ylim=ylim,"")
      }

    if(do.points) points(x$x,x$y)

    for (i in 1:n)
      {
        if(x$node[i])
          # Triangle i has positive area.
          # Connect circumcircle center of triangle i with neighbours:
          {
            # Find neighbour triangles
            tns<-sort(c(x$n1[i],x$n2[i],x$n3[i]))
            for(j in 1:3)
              {
                # Connect (if triangle exists and has positive area).
                if(tns[j]>0)
                  {
                  # simple node
                    if(x$node[tns[j]])
                      lines(c(x$x[i],x$x[tns[j]]),
                            c(x$y[i],x$y[tns[j]]),...)
                  }
                else if(tns[j]<0){
                  # dummy node
                  lines(c(x$x[i],x$dummy.x[-tns[j]]),
                        c(x$y[i],x$dummy.y[-tns[j]]),
                        lty="dashed",...) }
              }
          }
      }
    if(!add)
      title(main = main, sub =sub)
  }
