\name{tri.mesh}
\alias{tri.mesh}
\title{
  Create a delaunay triangulation
}
\usage{
tri.mesh(x, y = NULL, duplicate = "error",
         jitter = 10^-12, jitter.iter = 6, jitter.random = FALSE)
}
\arguments{
  \item{x}{vector containing x coordinates of the data. If \code{y} is missing
   \code{x} should contain two elements \code{$x} and \code{$y}.
  }
  \item{y}{vector containing y coordinates of the data.
  }
  \item{duplicate}{flag indicating how to handle duplicate elements.
   Possible values are: \code{"error"} -- default, \code{"strip"} -- remove all duplicate points, \code{"remove"} -- leave one point of duplicate points.
 }
   \item{jitter}{Jitter of amount of \code{diff(range(XX))*jitter} (XX=x
    or y) will be added to coordinates if collinear points are
    detected. Afterwards interpolation will be tried once again.

    Note that the jitter is not generated randomly unless
    \code{jitter.random} is set to \code{TRUE}. This ensures
    reproducable results. \code{\link[akima]{interp}} of package
    \code{akima} uses the same jitter mechanism. That means you can
    plot the triangulation on top of the interpolation and see the
    same triangulation as used for interpolation, see examples for
    \code{\link[akima]{interp}}.
  }
  \item{jitter.iter}{number of iterations to retry with jitter, amount
    will be increased in each iteration by \code{iter^1.5}}
  \item{jitter.random}{logical, see \code{jitter}, defaults to
    \code{FALSE}
  }
}
\value{
An object of class \code{"tri"}
}
\description{

  This subroutine creates a Delaunay triangulation of a
  set of N arbitrarily distributed points in the plane referred to as
  nodes.
  The Delaunay triangulation is defined
  as a set of triangles with the following five properties:

 1)  The triangle vertices are nodes.

 2)  No triangle contains a node other than its vertices.

 3)  The interiors of the triangles are pairwise disjoint.

 4)  The union of triangles is the convex hull of the set
       of nodes (the smallest convex set which contains
       the nodes).

 5)  The interior of the circumcircle of each triangle
       contains no node.

The first four properties define a triangulation, and the
last property results in a triangulation which is as close
as possible to equiangular in a certain sense and which is
uniquely defined unless four or more nodes lie on a common
circle.  This property makes the triangulation well-suited
for solving closest point problems and for triangle-based
interpolation.

  The triangulation can be generalized to a constrained
Delaunay triangulation by a call to \code{add.constraint}.
This allows for user-specified boundaries defining a nonconvex
and/or multiply connected region.

  The operation count for constructing the triangulation
is close to O(N) if the nodes are presorted on X or Y components.
Also, since the algorithm proceeds by adding
nodes incrementally, the triangulation may be updated with
the addition (or deletion) of a node very efficiently.
The adjacency information representing the triangulation
is stored as a linked list requiring approximately 13N
storage locations.
}
\references{
  R. J. Renka (1996). Algorithm 751: TRIPACK: a constrained
  two-dimensional Delaunay triangulation package.
  ACM Transactions on Mathematical Software.
  \bold{22}, 1-8.
}
\seealso{
  \code{\link{tri}}, \code{\link{print.tri}}, \code{\link{plot.tri}},
  \code{\link{summary.tri}}, \code{\link{triangles}},
  \code{\link{convex.hull}}, \code{\link{neighbours}},
  \code{\link{add.constraint}}.
}
\examples{
data(tritest)
tritest.tr<-tri.mesh(tritest$x,tritest$y)
tritest.tr
}
\keyword{spatial}
