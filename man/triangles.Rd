\name{triangles}
\title{Extract a list of triangles from a triangulation object}
\author{A. Gebhardt}
\usage{triangles(tri.obj)
}
\alias{triangles}
\arguments{
  \item{tri.obj}{object of class \code{"tri"}}
}
\description{This function extracts a triangulation data structure
  from an triangulation object created by \code{tri.mesh}.

    The vertices in the returned matrix (let's denote it  with
    \code{retval}) are ordered
    counterclockwise with the first vertex taken
    to be the one with smallest index.  Thus,
    \code{retval[i,"node2"]} and \code{retval[i,"node3"]} are larger
    than
    \code{retval[i,"node3"]} and index adjacent neighbors of
    node \code{retval[i,"node1"]}. The columns \code{trx} and
    \code{arcx}, x=1,2,3 index the triangle and arc,
    respectively, which are opposite (not shared
    by) node \code{nodex}, with \code{trix}= 0 if
    \code{arcx} indexes a boundary arc.  Vertex
    indexes range from 1 to N, triangle indexes
    from 0 to NT, and, if included, arc indexes
    from 1 to NA = NT+N-1.  The triangles are 
    ordered on first (smallest) vertex indexes,
    except that the sets of constraint triangles
    (triangles contained in the closure of a constraint
    region) follow the non-constraint
    triangles.
}
\value{A matrix with columns \code{node1},\code{node2},\code{node3},
  representing the vertex nodal indexes,
  \code{tr1},\code{tr2},\code{tr3}, representing neighboring triangle
  indexes and \code{arc1},\code{arc2},\code{arc3} reresenting arc indexes.

  Each row represents one triangle. 
}
\references{
  R. J. Renka (1996). Algorithm 751: TRIPACK: a constrained
  two-dimensional Delaunay triangulation package.
  ACM Transactions on Mathematical Software.
  \bold{22}, 1-8.
}
\seealso{\code{\link{tri}}, \code{\link{print.tri}}, \code{\link{plot.tri}}, \code{\link{summary.tri}}, \code{\link{triangles}}
}
\examples{
# use a slighlty modified version of data(tritest)
data(tritest2)
tritest2.tr<-tri.mesh(tritest2$x,tritest2$y)
triangles(tritest2.tr)
}
\keyword{spatial}
