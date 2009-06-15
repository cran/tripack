\name{plot.tri}
\title{Plot a triangulation object}
\author{A. Gebhardt}
\usage{\method{plot}{tri}(x, add=FALSE,xlim=range(x$x),ylim=range(x$y), do.points=TRUE, do.labels = FALSE, isometric=FALSE,...)
}
\alias{plot.tri}
\arguments{
  \item{x}{object of class \code{"tri"}}
  \item{add}{logical, if \code{TRUE}, add to a current plot.}
  \item{do.points}{logical, indicates if points should be plotted.}
  \item{do.labels}{logical, indicates if points should be labelled}
  \item{xlim,ylim}{x/y ranges for plot}
  \item{isometric}{generate an isometric plot (default \code{FALSE})}
  \item{...}{additional plot parameters}
}
\description{plots the triangulation \code{"x"}
}
\value{None
}
\references{
   R. J. Renka (1996). Algorithm 751: TRIPACK: a constrained
  two-dimensional Delaunay triangulation package.
  ACM Transactions on Mathematical Software.
  \bold{22}, 1-8. 
}
\seealso{
  \code{\link{tri}}, \code{\link{print.tri}}, 
  \code{\link{summary.tri}}
}
\examples{
# random points
plot(tri.mesh(rpois(100,lambda=20),rpois(100,lambda=20),duplicate="remove"))
# use a part of the quakes data set:
data(quakes)
quakes.part<-quakes[(quakes[,1]<=-10.78 & quakes[,1]>=-19.4 &
                     quakes[,2]<=182.29 & quakes[,2]>=165.77),]
quakes.tri<-tri.mesh(quakes.part$lon, quakes.part$lat, duplicate="remove")
plot(quakes.tri)
# use the whole quakes data set 
# (will not work with standard memory settings, hence commented out)
#plot(tri.mesh(quakes$lon, quakes$lat, duplicate="remove"), do.points=F) 
}
\keyword{spatial}
