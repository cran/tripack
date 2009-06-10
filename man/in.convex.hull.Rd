\name{in.convex.hull}
\title{Determines if points are in the convex hull of a triangulation object}
\author{A. Gebhardt}
\usage{in.convex.hull(tri.obj, x, y)
}
\alias{in.convex.hull}
\arguments{
  \item{tri.obj}{object of class \code{"tri"}}
  \item{x}{vector of x-coordinates of points to locate}
  \item{y}{vector of y-coordinates of points to locate}
}
\description{
Given a triangulation \code{tri.obj} of \eqn{n} points in the plane, this
subroutine returns a logical vector indicating if the points
\eqn{(x_i,y_i)}
are contained within the convex hull of \code{tri.obj}.

}
\value{
Logical vector.
}
\references{
  R. J. Renka (1996). Algorithm 751: TRIPACK: a constrained
  two-dimensional Delaunay triangulation package.
  ACM Transactions on Mathematical Software.
  \bold{22}, 1-8.
}
\seealso{
  \code{\link{tri}}, \code{\link{print.tri}}, \code{\link{plot.tri}},
  \code{\link{summary.tri}}, \code{\link{triangles}},
  \code{\link{add.constraint}}, \code{\link{convex.hull}}.
}
\examples{
# example from TRIPACK:
data(tritest)
tr<-tri.mesh(tritest$x,tritest$y)
in.convex.hull(tr,0.5,0.5)
in.convex.hull(tr,c(0.5,-1,1),c(0.5,1,1))
# use a part of the quakes data set:
data(quakes)
quakes.part<-quakes[(quakes[,1]<=-10.78 & quakes[,1]>=-19.4 &
                     quakes[,2]<=182.29 & quakes[,2]>=165.77),]
q.tri<-tri.mesh(quakes.part$lon, quakes.part$lat, duplicate="remove")
in.convex.hull(q.tri,quakes$lon[990:1000],quakes$lat[990:1000])
}
\keyword{spatial}
