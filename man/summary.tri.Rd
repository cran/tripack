\name{summary.tri}
\title{Return a summary of a triangulation object}
\author{A. Gebhardt}
\usage{summary(tri.obj)
}
\alias{summary.tri}
\arguments{
  \item{tri.obj}{object of class \code{"tri"}}
}
\description{
  Returns some information (number of nodes, triangles, arcs, boundary
  nodes and constraints) about \code{tri.obj}.
}
\value{An objekt of class \code{"summary.tri"}, to be printed by
  \code{\link{print.summary.tri}}.

  It contains the number of nodes (\code{n}), of arcs (\code{na}), of
  boundary nodes (\code{nb}), of triangles (\code{nt}) and constraints
  (\code{nc}).   
}
\references{
    R. J. Renka (1996). Algorithm 751: TRIPACK: a constrained
  two-dimensional {Delaunay} triangulation package.
  ACM Transactions on Mathematical Software.
  \bold{22}, 1-8.
}
\seealso{
  \code{\link{tri}}, \code{\link{print.tri}}, \code{\link{plot.tri}},
  \code{\link{print.summary.tri}}.
}
\keyword{spatial}
