\name{makezlist}
\alias{makezlist}
\title{Make a Zernike polynomial list}
\description{
  Makes a complete list of Zernike polynomial indices \code{(n, m, t)}
  suitable for subsequent calls to \code{\link{fillzm}}, \code{\link{pupil}}, etc.
}
\usage{
makezlist(minorder = 2, maxorder = 12)
}
\arguments{
  \item{minorder}{minimum radial polynomial order (must be even)}
  \item{maxorder}{maximum radial polynomial order (must be even)}
}
\value{
  A list with the following components:
  \item{n}{Radial polynomial order}
  \item{m}{Azimuthal order}
  \item{t}{character for trig function: one of c("n", "c", "s")}
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  This is a popular one dimensional ordering of the indices \code{(n,m)}
  of Zernike polynomials.

  The returned list is used in the functions \code{\link{fillzm}},
  \code{\link{pupil}}, \code{\link{zmult}}, etc.
}
\keyword{Mathematics}
