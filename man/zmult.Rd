\name{zmult}
\alias{zmult}
\title{Zernike coefficient multipliers}
\description{
  Conversion factors between normalized and conventional Zernike polynomials.
}
\usage{
zmult(zlist = zlist.qf)
}
\arguments{
  \item{zlist}{A list with named components \code{n, m, t}}
}
\details{
  The list indicates the Zernike polynomial orders to use, in the form
  returned by \code{\link{makezlist}} and {\code{\link{zlist.qf}}}.
}
\value{
  A vector the same length as the components of \code{zlist}.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  This function is not actually used in the current version of Rfringe, but is
  included for user convenience.
}
\seealso{
  \code{\link{Zernike}},
  \code{\link{makezlist}},
  \code{\link{zlist.qf}}.
}
\examples{
zcoef <- rnorm(length(zlist.qf$n))
zcoef		# a vector of normalized Zernike coefficients
zcoef*zmult()   # Coefficients in conventional representation
sqrt(crossprod(zcoef)) # This is the RMS error of the wavefront 
			# constructed from these Zernikes
}
\keyword{Mathematics}
