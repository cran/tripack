\name{tri.find}
\title{Locate a point in a triangulation}
\author{A. Gebhardt}
\usage{
tri.find(tri.obj,x,y)
}
\alias{tri.find}
\arguments{
  \item{tri.obj}{an triangulation object}
    \item{x}{x-coordinate of the point}
    \item{y}{y-coordinate of the point}
}
\description{
  This subroutine locates a point P=(\code{x},\code{y}) relative to a triangulation
  created by \code{tri.mesh}.  If P is
  contained in a triangle, the three vertex indexes are
  returned.  Otherwise, the indexes of the rightmost and
  leftmost visible boundary nodes are returned.
}
\value{
  A list with elements \code{i1},\code{i2},\code{i3} containing
  nodal indexes, in counterclockwise order,
  of the vertices of a triangle containing
  P=(\code{x},\code{y}), or, if P is not contained in the convex
  hull of the nodes, \code{i1} indexes the
  rightmost visible boundary node, \code{i2} indexes
  the leftmost visible boundary node,
  and \code{i3} = 0.  Rightmost and leftmost are
  defined from the perspective of P, and a
  pair of points are visible from each
  other if and only if the line segment
  joining them intersects no triangulation
  arc.  If P and all of the nodes lie on a
  common line, then \code{i1}=\code{i2}=\code{i3} = 0 on
  output.
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
  \code{\link{convex.hull}}
}
\examples{
data(tritest)
tritest.tr<-tri.mesh(tritest$x,tritest$y)
plot(tritest.tr)
pnt<-list(x=0.3,y=0.4)
triangle.with.pnt<-tri.find(tritest.tr,pnt$x,pnt$y)
attach(triangle.with.pnt)
lines(tritest$x[c(i1,i2,i3,i1)],tritest$y[c(i1,i2,i3,i1)],col="red")
points(pnt$x,pnt$y)
}
\keyword{spatial}
