\name{plot.voronoi}
\title{Plot a xect}
\author{A. Gebhardt}
\usage{\method{plot}{voronoi}(x, add=FALSE, xlim=c(min(x$tri$x)-0.1*diff(range(x$tri$x)), max(x$tri$x) + 0.1 *diff(range(x$tri$x))),   ylim=c(min(x$tri$y)-0.1*diff(range(x$tri$y)),max(x$tri$y)+0.1*diff(range(x$tri$y))),all=FALSE,do.points=TRUE,main="Voronoi mosaic",sub=deparse(substitute(x)), ...)
}
\alias{plot.voronoi}
\arguments{
  \item{x}{object of class \code{"voronoi"}}
  \item{add}{logical, if \code{TRUE}, add to a current plot.}
  \item{xlim}{x plot ranges, by default modified to hide dummy 
        points outside of the plot}
  \item{ylim}{y plot ranges, by default modified to hide dummy 
        points outside of the plot}
  \item{all}{show all (including dummy points in the plot}
  \item{do.points}{logical, indicates if points should be plotted.}
  \item{main}{plot title}
  \item{sub}{plot subtitle}
  \item{...}{additional plot parameters}
}
\description{Plots the mosaic \code{"x"}.

  Dashed lines are used for outer tiles of the mosaic.
}
\value{None
}
\references{
   R. J. Renka (1996). Algorithm 751: TRIPACK: a constrained
  two-dimensional {Delaunay} triangulation package.
  ACM Transactions on Mathematical Software.
  \bold{22}, 1-8. 
}
\seealso{
  \code{\link{voronoi}}, \code{\link{print.voronoi}}, 
  \code{\link{summary.voronoi}}
}
\examples{
# generate a random mosaic
plot(voronoi.mosaic(runif(100),runif(100),duplicate="remove"))
# use a part of the quakes data set:
data(quakes)
quakes.part<-quakes[(quakes[,1]<=-17 & quakes[,1]>=-19.0 &
                     quakes[,2]<=182.0 & quakes[,2]>=180.0),]
quakes.vm<-voronoi.mosaic(quakes.part$lon, quakes.part$lat, duplicate="remove")
plot(quakes.vm)
# use the whole quakes data set 
# (will not work with standard memory settings, hence commented out here)
#plot(voronoi.mosaic(quakes$lon, quakes$lat, duplicate="remove"))
}
\keyword{spatial}
