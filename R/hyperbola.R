hyperbola <- function(x1,y1,x2,y2,w1,w2,tl,tu,...){


  transrot <- function(x,y,x0,y0,rho){
    n <- length(x)
    Atrans <- matrix(c(1,0,x0,
                       0,1,y0,
                       0,0, 1),
                     nrow=3,ncol=3,byrow=TRUE)

    Arot <- matrix(c(cos(rho) ,sin(rho),0,
                     -sin(rho),cos(rho),0,
                             0,       0,1),
                   nrow=3,ncol=3,byrow=TRUE)

    ret <-  Atrans %*% Arot %*% t(cbind(x,y,rep(1,n)))
           
    list(x=ret[1,],y=ret[2,])
  }
  
  t <- seq(tl,tu,length=50)
  cx <- (x1+x2)/2
  cy <- (y1+y2)/2
  phi <- -atan2(y2-y1,x2-x1)

  a <- (w2-w1)/2
  d <- sqrt((x2-x1)^2+(y2-y1)^2)
  b <- sqrt((d/2)^2-a^2)

  
  hx <- -a*cosh(t)
  hy <- b*sinh(t)


  h <- transrot(hx,hy,cx,cy,phi)

  lines(h$x,h$y,...)
#  text(h$x[15],h$y[15],tl)
#  text(h$x[35],h$y[35],tu)
  
}
