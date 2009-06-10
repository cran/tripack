\name{voronoi.polygons}
\alias{voronoi.polygons}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ extract polygons from a voronoi mosaic }
\description{
  This functions extracts polygons from a \code{voronoi.mosaic} object.
}
\usage{
voronoi.polygons(voronoi.obj)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{voronoi.obj}{ object of class \code{voronoi.mosaic} }
}
\value{
  Returns an object of class \code{voronoi.polygons} with unamed list
  elements for each polygon. These list
  elements are matrices with columns \code{x} and \code{y}.
}
\author{ Denis White }
\seealso{ \code{\link{plot.voronoi.polygons}},\code{\link{voronoi.mosaic}}}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	     or do  help(data=index)  for the standard data sets.

data(tritest)
tritest.vm <- voronoi.mosaic(tritest$x,tritest$y)
tritest.vp <- voronoi.polygons(tritest.vm)
tritest.vp
}
\keyword{ spatial }
