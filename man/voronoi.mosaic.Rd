\name{voronoi.mosaic}
\title{Create a Voronoi mosaic}
\author{A. Gebhardt}
\usage{voronoi.mosaic(x,y=NULL,duplicate="error")
}
\alias{voronoi.mosaic}
\arguments{
  \item{x}{vector containing x coordinates of the data. If \code{y} is
    missing    \code{x} should contain two elements \code{$x} and \code{$y}.
  }
  \item{y}{vector containing y coordinates of the data.
  }
  \item{duplicate}{flag indicating how to handle duplicate elements.
   Possible values are: \code{"error"} -- default, \code{"strip"} --
   remove all duplicate points, \code{"remove"} -- leave one point of
   duplicate points. 
 }
}
\description{
  This function creates a Voronoi mosaic.

  It creates first a Delaunay triangulation, determines the circumcircle 
  centers of its triangles, and connects these points according to the
  neighbourhood relations between the triangles.
}
\value{
  An object of class \code{\link{voronoi}}.
}
\seealso{  \code{\link{voronoi}},\code{\link{voronoi.mosaic}},      \code{\link{print.voronoi}}, \code{\link{plot.voronoi}}
}
\examples{
# example from TRIPACK:
data(tritest)
tritest.vm<-voronoi.mosaic(tritest$x,tritest$y)
tritest.vm
# use a part of the quakes data set:
quakes.part<-quakes[(quakes[,1]<=-17 & quakes[,1]>=-19.0 &
                     quakes[,2]<=182.0 & quakes[,2]>=180.0),]
quakes.vm<-voronoi.mosaic(quakes.part$lon, quakes.part$lat, duplicate="remove")
quakes.vm
}
\keyword{}
