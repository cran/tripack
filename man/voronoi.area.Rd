\name{voronoi.area}
\title{Calculate area of Voronoi polygons}
\author{S. J. Eglen}
\usage{voronoi.area(voronoi.obj)
}
\alias{voronoi.area}
\arguments{
  \item{voronoi.obj}{object of class \code{"voronoi"}}
}
\description{Computes the area of each Voronoi polygon.
For some sites at the edge of the region, the Voronoi polygon is not
bounded, and so the area of those sites cannot be calculated, and hence
will be NA.
}
\value{A vector of polygon areas.}
\seealso{
  \code{\link{voronoi}},
}
\keyword{spatial}

\examples{
data(tritest)
tritest.vm <- voronoi.mosaic(tritest$x,tritest$y)
tritest.vm.areas <- voronoi.area(tritest.vm)
plot(tritest.vm)
text(tritest$x, tritest$y, tritest.vm.areas)
}

