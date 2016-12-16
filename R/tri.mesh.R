tri.mesh <- function(x,y=NULL,duplicate="error",
                     jitter = 10^-12, jitter.iter = 6, jitter.random = FALSE)
{
  if(is.null(x))
     stop("argument x missing.")
  if(is.null(y)){
    x1<-x$x
    y1<-x$y
    if (is.null(x1) || is.null(y1))
      stop("argument y missing and x contains no $x or $y component.")
  }
  else{
    x1<-x
    y1<-y
  }

  n <- length(x1)
  if(length(y1)!=n)
    stop("length of x and y differ.")
  # handle duplicate points:
  xy <- paste(x1, y1, sep =",")
  i <- match(xy, xy)
  if(duplicate!="error")
    {
      if(duplicate!="remove" & duplicate!="error" & duplicate!="strip"){
        stop("possible values for \'duplicate\' are \"error\", \"strip\" and \"remove\"")
      }
      else{
        if(duplicate=="remove")
          ord <- !duplicated(xy)
        if(duplicate=="strip")
          ord <- (hist(i,plot=FALSE,freq=TRUE,breaks=seq(0.5,max(i)+0.5,1))$counts==1)
        x1 <- x1[ord]
        y1 <- y1[ord]
        n <- length(x1)
      }
    }
  else
    if(any(duplicated(xy)))
      stop("duplicate data points")

  ans<-.Fortran("trmesh",
                as.integer(n),
                x=as.double(x1),
                y=as.double(y1),
                tlist=integer(6*n-12),
                tlptr=integer(6*n-12),
                tlend=integer(n),
                tlnew=as.integer(0),
                tnear=integer(n),
                tnext=integer(n),
                tdist=double(n),
                ier=as.integer(0),
                PACKAGE = "tripack")
  if(ans$ier==0)
  {
      tri.obj<-list(n=n,x=x1,y=y1,tlist=ans$tlist,tlptr=ans$tlptr,
                    tlend=ans$tlend,tlnew=ans$tlnew,
                    nc=0,lc=0,call=match.call())
  } else if(ans$ier==-2) {
      ## first 3 nodes collinear, retry with added jitter to avoid colinearities!
      jitter.trials <- 1
      success <- FALSE
      while(jitter.trials<jitter.iter & !success){
          ## dont use random jitter for reproducabilty by default,
          ## same jitter as in akima interp() function:
          if(jitter.random){
              ## the randomness is only contained in a random shift
              ## of a regular -1,+/-0,+1,-1,+/-0,+1,... scheme
              ## determining the fixed amounts of jitter to be added
              ## to the x and y values separately.
              ## Really random jitter like this
              ## xj <- x+rnorm(n,0,diff(range(x1))*jitter*jitter.trials^1.5)
              ## yj <- y+rnorm(n,0,diff(range(y1))*jitter*jitter.trials^1.5)
              ## still triggers spurious errors in the triangulation.
              ## maybe reordering would also suffice:
              ## j <- sample(1:length(x), length(x))
              ## xj <- x[j]
              ## yj <- y[j]
              ## but it needs afterwards back ordering!
              j <- list()
              j[[1]] <- rep(c(-1,0,1),length.out=length(x1))
              j[[2]] <- rep(c(0,1,-1),length.out=length(x1))
              j[[3]] <- rep(c(1,-1,0),length.out=length(x1))
              jx <- sample(1:3,1)
              jy <- sample(1:3,1)
              xj <- x1+j[[jx]]*diff(range(x1))*jitter*jitter.trials^1.5
              yj <- y1+j[[jy]]*diff(range(y1))*jitter*jitter.trials^1.5
          } else {
              xj <- x1+rep(c(-1,0,1),length.out=length(x1))*diff(range(x1))*jitter*jitter.trials^1.5
              yj <- y1+rep(c(0,1,-1),length.out=length(y1))*diff(range(y1))*jitter*jitter.trials^1.5
          }
          ans<-.Fortran("trmesh",
                        n=as.integer(n),
                        x=as.double(xj),
                        y=as.double(yj),
                        tlist=integer(6*n-12),
                        tlptr=integer(6*n-12),
                        tlend=integer(n),
                        tlnew=as.integer(0),
                        tnear=integer(n),
                        tnext=integer(n),
                        tdist=double(n),
                        ier=as.integer(0),
                        PACKAGE = "tripack")
          success <- (ans$ier==0)
          if(success)
              warning("success: collinearities reduced through jitter")

          jitter.trials <- jitter.trials+1
      }
      if(ans$ier==0){
          warning("dataset started with 3 colinear points, jitter added!")
          tri.obj<-list(n=n,x=x1,y=y1,tlist=ans$tlist,tlptr=ans$tlptr,
                        tlend=ans$tlend,tlnew=ans$tlnew,
                        nc=0,lc=0,call=match.call())
      } else
          stop(paste("error in trmesh, code: ",ans$ier,", collinear points."))

  } else
      stop(paste("error in trmesh, code: ",ans$ier))

  class(tri.obj)<-"tri"
  invisible(tri.obj)
}
