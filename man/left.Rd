\name{left}
\Rdversion{1.1}
\alias{left}
\alias{is.left.of}
\title{
  Determines whether given points are left of a directed edge.
}
\description{
  This function returns a logical vector indicating which elements of
  the given points P0 are left of the directed edge P1->P2.
}
\usage{
left(x0, y0, x1, y1, x2, y2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x0}{
    Numeric vector, 'x' coordinates of points P0 to check
  }
  \item{y0}{
    Numeric vector, 'y' coordinates of points P0 to check, same length as 'x'.
  }
  \item{x1}{
    'x' coordinate of point P1
  }
  \item{y1}{
    'y' coordinate of point P1
  }
  \item{x2}{
    'x' coordinate of point P2
  }
  \item{y2}{
    'y' coordinate of point P2
  }
}

\value{
  Logical vector.
}
\author{
  A. Gebhardt
}
\note{
  This is an interface to the Fortran function VLEFT, wich is modeled
  after TRIPACKs LEFT function but accepts more than one point P0.
}

\seealso{
  \code{\link{in.convex.hull}}
}
\examples{

left(c(0,0,1,1),c(0,1,0,1),0,0,1,1)
   
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }

