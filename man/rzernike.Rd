\name{rzernike}
\alias{rzernike}
\alias{odd}
\title{Radial Zernike Polynomial}
\description{Zernike's Radial Polynomials}
\usage{
rzernike(rho, n, m)
odd(n,m)
}
\arguments{
  \item{rho}{normalized radius, \eqn{0 <= rho <= 1}}
  \item{n}{radial polynomial order}
  \item{m}{azimuthal order}
}
\details{
  The arguments \code{n} and \code{m} must be relatively even.
}
\value{
  The value of the radial Zernike polynomial of order \code{(n,m)} at
  normalized radius \code{rho}. If \code{rho} is a vector, matrix, or higher order array
  the returned value is a vector or array of the same dimension.
  
  The utility function \code{odd(n,m)} returns \code{FALSE} iff \code{(n-m)\%\%2 == 0}.
}
\references{
  \url{http://mathworld.wolfram.com/ZernikePolynomial.html}
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  This routine implements the recurrence relation given in Equation 12 of
  \url{http://mathworld.wolfram.com/ZernikePolynomial.html}.

  In general you should call the higher level function \code{\link{Zernike}} instead
  of this.
}
\seealso{
  \code{\link{Zernike}}.
}
\keyword{Mathematics}
