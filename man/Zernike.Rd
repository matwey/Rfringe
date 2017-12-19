\name{Zernike}
\alias{Zernike}
\title{Zernike Polynomials}
\description{
  Routines for creating and manipulating Zernike polynomials.
}
\usage{
Zernike(rho, theta, n, m, t)
}
\arguments{
  \item{rho}{normalized radius, \eqn{0 <= rho <= 1}}
  \item{theta}{angular coordinate}
  \item{n}{radial polynomial order}
  \item{m}{azimuthal order}
  \item{t}{character for trig function: one of c("n", "c", "s")}
}
\details{
  The arguments \code{n} and \code{m} must be relatively even.
}
\value{
  The value of the Zernike polynomial of order \code{(n, m)}
  at polar coordinates \code{(rho, theta)}. The arguments
  \code{rho} and \code{theta} may be vectors, matrices, or higher order arrays, in which case
  the returned value is a vector or array of the same dimension.
}
\references{
  Born, M. and Wolf, E. 1999, \emph{Principles of Optics, 7th Edition},
  Cambridge University Press, chapter 9 and appendix VII.

  Noll, R.J. 1976, \bold{Zernike polynomials and atmospheric turbulence}, \emph{J. Opt. Soc. Am.},
  Vol. 66, No. 3, p. 207.

  \url{http://wyant.opt-sci.arizona.edu/zernikes/zernikes.htm}

  \url{http://mathworld.wolfram.com/ZernikePolynomial.html}
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  This function returns Zernikes scaled such that they form an orthonormal basis set
  for the space of functions defined on the unit circle. Note that this is not the most
  commonly used definition (as given e.g. in \cite{Born and Wolf}). The definition I use is
  often associated with \cite{Noll (1976)}.

  The otherwise unused function \code{\link{zmult}} can be used to convert between normalized
  and conventionally defined vectors of Zernike coefficients.
}
\seealso{
  \code{\link{rzernike}},
  \code{\link{makezlist}},
  \code{\link{zlist.qf}},
  \code{\link{zmult}},
  \code{\link{fillzm}},
  \code{\link{pupil}},
  \code{\link{pupilrms}},
  \code{\link{pupilpv}},
  \code{\link{strehlratio}}.
}
\examples{
Zernike(1, 0, 4, 0, "n")	# == sqrt(5)

# A slightly more complex example

rho <- seq(0, 1, length = 101)
theta <- rep(0, 101)

plot(rho, Zernike(rho, theta, 6, 0, "n"), type="l", 
  ylim=c(-3.5,3.5), main="Some 6th order Zernike Polynomials")
lines(rho, Zernike(rho, theta, 5, 1, "c"), lty=2)
lines(rho, Zernike(rho, theta, 4, 2, "c"), lty=3)
lines(rho, Zernike(rho, theta, 3, 3, "c"), lty=4)
}
\keyword{Mathematics}
