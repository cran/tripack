\name{add.constraint}
\title{
 Add a constraint to an triangulaion object
}
\alias{add.constraint}
\usage{
add.constraint(tri.obj,cstx,csty,reverse=F)
}
\arguments{
  \item{tri.obj}{object of class \code{"tri"}}
  \item{cstx}{vector containing x coordinates of the constraint curve.}
  \item{csty}{vector containing y coordinates of the constraint curve.}
  \item{reverse}{if \code{TRUE} the orientation of the constraint curve
    is reversed.}
}
\value{
An new object of class \code{"tri"}.}
#\notes{
#}
\description{
  This subroutine provides for creation of a constrained
Delaunay triangulation which, in some sense, covers an
arbitrary connected region R rather than the convex hull
of the nodes.  This is achieved simply by forcing the
presence of certain adjacencies (triangulation arcs) corresponding
to constraint curves.  The union of triangles
coincides with the convex hull of the nodes, but triangles
in R can be distinguished from those outside of R.  The
only modification required to generalize the definition of
the Delaunay triangulation is replacement of property 5
(refer to \code{\link{tri.mesh}} by the following:

5')  If a node is contained in the interior of the circumcircle
of a triangle, then every interior point
of the triangle is separated from the node by a
constraint arc.

In order to be explicit, we make the following definitions.
A constraint region is the open interior of a
simple closed positively oriented polygonal curve defined
by an ordered sequence of three or more distinct nodes
(constraint nodes) P(1),P(2),...,P(K), such that P(I) is
adjacent to P(I+1) for I = 1,...,K with P(K+1) = P(1).
Thus, the constraint region is on the left (and may have
nonfinite area) as the sequence of constraint nodes is
traversed in the specified order.  The constraint regions
must not contain nodes and must not overlap.  The region
R is the convex hull of the nodes with constraint regions
excluded.

  Note that the terms boundary node and boundary arc are
reserved for nodes and arcs on the boundary of the convex
hull of the nodes.

  The algorithm is as follows:  given a triangulation
which includes one or more sets of constraint nodes, the
corresponding adjacencies (constraint arcs) are forced to
be present (Fortran subroutine EDGE).  Any additional new arcs
required are chosen to be locally optimal (satisfy the
modified circumcircle property).
}
\references{
  R. J. Renka (1996). Algorithm 751: TRIPACK: a constrained
  two-dimensional {Delaunay} triangulation package.
  ACM Transactions on Mathematical Software.
  \bold{22}, 1-8.
}

\seealso{
  \code{\link{tri}}, \code{\link{print.tri}}, \code{\link{plot.tri}}, \code{\link{summary.tri}}, \code{\link{triangles}}, \code{\link{convex.hull}}.
}
\examples{
# we will use the simple test data from TRIPACK:
data(tritest)
tritest.tr<-tri.mesh(tritest)
opar<-par(mfrow=c(2,2))
plot(tritest.tr)
# include all points in a big triangle:
tritest.tr<-add.constraint(tritest.tr,c(-0.1,2,-0.1),
                                      c(-3,0.5,3),reverse=T)
# insert a small cube:
tritest.tr <- add.constraint(tritest.tr, c(0.4, 0.4,0.6, 0.6), 
                                         c(0.6, 0.4,0.4, 0.6), 
                                         reverse = F)
par(opar)
}
\keyword{}
