\name{tri.dellens}
\title{Compute the Delaunay segment lengths}
\author{S. J. Eglen}
\usage{
tri.dellens(voronoi.obj, exceptions = NULL, inverse = FALSE)
}
\alias{tri.dellens}
\arguments{
  \item{voronoi.obj}{object of class \code{"voronoi"}}
  \item{exceptions}{a numerical vector}
  \item{inverse}{Logical}
}
\description{
  Return a vector of Delaunay segment lengths for the voronoi object.
  The Delaunay triangles connected to sites contained in  \code{exceptions}
  vector are ignored (unless \code{inverse} is TRUE, when only those
  Delaunay triangles are accepted).

  The \code{exceptions} vector is provided so that sites at the border
  of a region can be removed, as these tend to bias the distribution of
  Delaunay segment lengths.  \code{exceptions} can be created by
  \code{\link{voronoi.findrejectsites}}.
}
\value{
  A vector of Delaunay segment lengths.
}
\seealso{
  \code{\link{voronoi.findrejectsites}}, \code{\link{voronoi.mosaic}},
}
\examples{
data(tritest)
tritest.vm <- voronoi.mosaic(tritest$x,tritest$y)

tritest.vm.rejects <- voronoi.findrejectsites(tritest.vm, 0,1, 0, 1)
trilens.all <- tri.dellens(tritest.vm)
trilens.acc <- tri.dellens(tritest.vm, tritest.vm.rejects)
trilens.rej <- tri.dellens(tritest.vm, tritest.vm.rejects, inverse=T)

par(mfrow=c(3,1))
dotchart(trilens.all, main="all Delaunay segment lengths")
dotchart(trilens.acc, main="excluding border sites")
dotchart(trilens.rej, main="only border sites")
}
\keyword{spatial}
