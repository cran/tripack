\name{identify.tri}
\title{Identify points in a triangulation plot}
\author{A. Gebhardt}
\usage{\method{identify}{tri}(x,...)
}
\alias{identify.tri}
\arguments{
  \item{x}{object of class \code{"tri"}}
  \item{...}{additional paramters for \code{identify}}
}
\description{Identify points in a plot of \code{"x"} with its
  coordinates. The plot of \code{"x"} must be generated with \code{plot.tri}.
}
\value{an integer vector containing the indexes of the identified points.
}
\seealso{
    \code{\link{tri}}, \code{\link{print.tri}}, \code{\link{plot.tri}}, \code{\link{summary.tri}}
}
\examples{
data(tritest)
tritest.tr<-tri.mesh(tritest$x,tritest$y)
plot(tritest.tr)
identify.tri(tritest.tr)
}
\keyword{spatial}
