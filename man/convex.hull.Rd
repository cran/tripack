\name{convex.hull}
\title{Return the convex hull of a triangulation object}
\author{A. Gebhardt}
\usage{convex.hull(tri.obj, plot.it=F, add=F,...)
}
\alias{convex.hull}
\arguments{
  \item{tri.obj}{object of class \code{"tri"}}
  \item{plot.it}{logical, if \code{TRUE} the convex hull of
    \code{tri.obj} will be plotted.}
  \item{add}{logical. if \code{TRUE} (and \code{plot.it=T}), add to
    a current plot.}
}
\description{
Given a triangulation \code{tri.obj} of \eqn{n} points in the plane, this
subroutine returns two vectors containing the coordinates
of the nodes on the boundary of the convex hull.

}
\value{
  \item{x}{x coordinates of boundary nodes.}
  \item{y}{y coordinates of boundary nodes.}
}
\references{
  R. J. Renka (1996). Algorithm 751: TRIPACK: a constrained
  two-dimensional {Delaunay} triangulation package.
  ACM Transactions on Mathematical Software.
  \bold{22}, 1-8.
}
\seealso{
  \code{\link{tri}}, \code{\link{print.tri}}, \code{\link{plot.tri}}, \code{\link{summary.tri}}, \code{\link{triangles}}, \code{\link{add.constraint}}.
}
\examples{
# rather simple example from TRIPACK:
tr<-tri.mesh(tritest$x,tritest$y)
convex.hull(tr,plot.it=T)
# random points:
rand.tr<-tri.mesh(runif(10),runif(10))
plot(rand.tr)
rand.ch<-convex.hull(rand.tr, plot.it=T, add=T, col="red")
# use a part of the quakes data set:
quakes.part<-quakes[(quakes[,1]<=-17 & quakes[,1]>=-19.0 &
                     quakes[,2]<=182.0 & quakes[,2]>=180.0),]
quakes.tri<-tri.mesh(quakes.part$lon, quakes.part$lat, duplicate="remove")
plot(quakes.tri)
convex.hull(quakes.tri, plot.it=T, add=T, col="red")
}
\keyword{}
