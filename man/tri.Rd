\name{tri}
\title{A triangulation object}
\author{A. Gebhardt}
\alias{tri}
\arguments{
  \item{n}{Number of nodes}
  \item{x}{x coordinates of the triangulation nodes}
  \item{y}{y coordinates of the triangulation nodes}}
  \item{tlist}{Set of nodal indexes which, along with \code{tlptr},
    \code{tlend}, and \code{tlnew}, define the triangulation as a
    set of \eqn{n} adjacency lists -- counterclockwise-ordered sequences of
    neighboring nodes such that the first and last neighbors of a
    boundary node are boundary nodes (the first neighbor of an interior
    node is arbitrary). In order to distinguish between
    interior and boundary nodes, the last neighbor of each boundary node
    is represented by the negative of its index.}
  \item{tlptr}{Set of pointers in one-to-one
    correspondence with the elements of \code{tlist}.
    \code{tlist[tlptr[i]]} indexes the node which follows
    \code{tlist[i]} in cyclical counterclockwise order
    (the first neighbor follows the last neighbor).}
  \item{tlend}{Set of pointers to adjacency lists. \code{tlend[k]}
    points to the last neighbor of node \eqn{k} for
    \eqn{k = 1,...,n}.  Thus, \code{tlist[\tlend[k]]}<0 if and
    only if \eqn{k} is a boundary node.}
  \item{tlnew}{ Pointer to the first empty location in \code{tlist}
    and \code{tlptr} (list length plus one).}
  \item{nc}{number of constraints}
  \item{lc}{starting indices of constraints in \code{x} and \code{y}}
  \item{call}{call, which generated this object}
}
\notes{
The elements \code{tlist}, \code{tlptr}, \code{tlend} and \code{tlnew}
are mainly intended for internal use in the appropriate Fortran
routines. 
}
\references{
    R. J. Renka (1996). Algorithm 751: TRIPACK: a constrained
  two-dimensional {Delaunay} triangulation package.
  ACM Transactions on Mathematical Software.
  \bold{22}, 1-8.
}
\seealso{\code{\link{tri.mesh}}, \code{\link{print.tri}}, \code{\link{plot.tri}}, \code{\link{summary.tri}}
}
\keyword{}
