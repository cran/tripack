\name{voronoi}
\title{Voronoi object}
\author{A. Gebhardt}
\alias{voronoi}
\arguments{
  \item{x,y}{x and y coordinates of nodes of the voronoi mosaic. Each node is a
    circumcircle center of some triangle from the Delaunay triangulation.}
  \item{node}{logical vector, indicating real nodes of the voronoi
  mosaic. These nodes are the centers of circumcircles of triangles with 
  positive area of the delaunay triangulation.

  If \code{node[i]}=\code{FALSE}, (\code{c[i]},\code{x[i]}) belongs to a 
  triangle with area 0.}
\item{n1,n2,n3}{indices of neighbour nodes. Negative indices indicate
  dummy points as neighbours.}
\item{tri}{triangulation object, see \code{\link{tri}}.}
\item{area}{area of triangle i. \code{area[i]}=-1 indicates a removed
  triangle with area 0 at the border of the triangulation.}
\item{ratio}{aspect ratio (inscribed radius/circumradius) of triangle
  i.}
\item{radius}{circumradius of triangle i.}
\item{dummy.x,dummy.y}{x and y coordinates of dummy points. They are
  used for plotting of unbounded tiles.}
}
\description{
  An \code{voronoi} object is created with \code{\link{voronoi.mosaic}
}}
\seealso{\code{\link{voronoi.mosaic}},\code{\link{plot.voronoi}}
}
\keyword{spatial}

