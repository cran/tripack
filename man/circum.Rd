\name{circum}
\Rdversion{1.1}
\alias{circum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Determine the circumcircle of a triangle
}
\description{
This function returns the circumcircle of a triangle.
}
\usage{
circum(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of three elements, giving the x coordinatres of the
    triangle nodes.
}
  \item{y}{Vector of three elements, giving the y coordinatres of the
    triangle nodes.
}
}
\details{
  This is an interface to the Fortran function CIRCUM found in TRIPACK.
}
\value{
  \item{ x }{ 'x' coordinate of center}
  \item{ y }{ 'y' coordinate of center}
  \item{ radius }{ circumcircle radius }
  \item{ signed.area }{ signed area of riangle (positive iff nodes are
    numbered counter clock wise) }
  \item{ aspect.ratio }{ ratio "radius of inscribed circle"/"radius of
    circumcircle", varies between 0 and 0.5

    0 means collinear points, 0.5 equilateral trangle.
    }
}
\references{
  R. J. Renka (1996). Algorithm 751: TRIPACK: a constrained
  two-dimensional Delaunay triangulation package.
  ACM Transactions on Mathematical Software.
  \bold{22}, 1-8.
}
\author{
  Fortran code: R. J. Renka, R code: A. Gebhardt
}
\note{
  This function is mainly intended to be used by \code{\link{circumcircle}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{circumcircle}}
}
\examples{
circum(c(0,1,0),c(0,0,1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }

