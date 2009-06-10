\name{summary.voronoi}
\title{Return a summary of a objectect}
\author{A. Gebhardt}
\usage{\method{summary}{voronoi}(object,...)
}
\alias{summary.voronoi}
\arguments{
  \item{object}{object of class \code{"voronoi"}}
  \item{...}{additional paramters for \code{summary}}
}
\description{
  Returns some information about \code{object}
}
\value{Object of class \code{"summary.voronoi"}.

  It contains the number of nodes (\code{nn}) and dummy nodes (\code{nd}).
}
\references{
    R. J. Renka (1996). Algorithm 751: TRIPACK: a constrained
  two-dimensional Delaunay triangulation package.
  ACM Transactions on Mathematical Software.
  \bold{22}, 1-8.
}
\seealso{
  \code{\link{voronoi}},\code{\link{voronoi.mosaic}},
  \code{\link{print.voronoi}}, \code{\link{plot.voronoi}},
  \code{\link{print.summary.voronoi}}. 
}
\keyword{spatial}
