\name{pupil}
\alias{pupil}
\title{Create a circular pupil and fill it with a wavefront}
\description{
  Creates a representation of a wavefront from a vector of Zernike polynomial values.
}
\usage{
pupil(size = 255, obstruct = 0, zcoef = NULL, zlist = zlist.qf, phi = 0, piston = 0)
}
\arguments{
  \item{size}{size of the returned matrix}
  \item{obstruct}{central obstruction \emph{fraction}}
  \item{zcoef}{Vector of Zernike coefficients}
  \item{zlist}{A list with named components \code{n, m, t} describing the contents of \code{zcoef}}
  \item{phi}{angular coordinate to rotate entire coordinate system}
  \item{piston}{Piston (constant) term to add to the wavefront}
}
\details{
  The coordinate system is rotated \emph{clockwise} by the angle \code{phi}
  specified in \emph{degrees}. This is done to present consistent displays
  of rotated wavefronts in \code{\link{Rfringe}}.
}
\value{
  A \code{size} by \code{size} matrix of wavefront values. \code{NA}'s are used to
  fill out the matrix outside the circular pupil and inside the obstruction.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  \emph{Most} high level \R graphics functions will handle \code{NA}'s as intended.

  This function can take a while if \code{zcoef} includes high order Zernikes. Decreasing
  the matrix size will help with speed, but may provide too low resolution for good graphical
  representations. No attempt is made to \dQuote{anti-alias} the edges of the pupil.
}
\seealso{
  \code{\link{Zernike}},
  \code{\link{makezlist}},
  \code{\link{zlist.qf}},
  \code{\link{fillzm}},
  \code{\link{pupilrms}},
  \code{\link{pupilpv}}.
}
\examples{
# A fairly typical use of this function:

wf <- pupil(zcoef=rnorm(length(zlist.qf$n)))
image(wf, col=topo.colors(256), asp=1)
contour(wf, add=TRUE)
}
\keyword{Mathematics}
