\name{neighbours}
\title{List of neighbours from a triangulation object}
\author{A. Gebhardt}
\usage{neighbours(tri.obj)
}
\alias{neighbours}
\arguments{
  \item{tri.obj}{object of class \code{"tri"}}
}
\description{Extract a list of neighbours from a triangulation object
}
\value{nested list of neighbours per point
}
\references{
  R. J. Renka (1996). Algorithm 751: TRIPACK: a constrained
  two-dimensional {Delaunay} triangulation package.
  ACM Transactions on Mathematical Software.
  \bold{22}, 1-8.
}
\seealso{\code{\link{tri}}, \code{\link{print.tri}}, \code{\link{plot.tri}}, \code{\link{summary.tri}}, \code{\link{triangles}}
}
\examples{
data(tritest)
tritest.tr<-tri.mesh(tritest$x,tritest$y)
tritest.nb<-neighbours(tritest.tr)
}
\keyword{}
