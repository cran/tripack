\name{plot.voronoi.polygons}
\alias{plot.voronoi.polygons}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{plots an voronoi.polygons object}
\description{
  plots an \code{voronoi.polygons} object
}
\usage{
plot.voronoi.polygons(x, which, color=TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ object of class \code{voronoi.polygons}  }
  \item{which}{ index vector selecting which polygons to plot }
  \item{color}{ logical, determines if plot should be colored, default: \code{TRUE} }
  \item{\dots}{ additional plot arguments }
}
\author{ A. Gebhardt}

\seealso{ \code{\link{voronoi.polygons}}}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	     or do  help(data=index)  for the standard data sets.
data(tritest)
tritest.vm <- voronoi.mosaic(tritest$x,tritest$y)
tritest.vp <- voronoi.polygons(tritest.vm)
plot(tritest.vp)
plot(tritest.vp,which=c(1,3,5))
}
\keyword{ spatial }% at least one, from doc/KEYWORDS
