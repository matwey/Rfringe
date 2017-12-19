\name{synth.interferogram}
\alias{synth.interferogram}
\title{Synthetic Interferogram}
\description{
  Computes and displays a synthetic interferogram for a wavefront constructed
  from a vector of Zernike coefficients.
}
\usage{
synth.interferogram(zcoef, zlist = zlist.qf, phi = 0, 
  size = 255, obstruct = 0, iname = "")
}
\arguments{
  \item{zcoef}{Vector of Zernike coefficients, with piston term as the first element}
  \item{zlist}{A list with named components \code{n, m, t} describing the contents of \code{zcoef}}
  \item{phi}{angular coordinate to rotate entire coordinate system}
  \item{size}{Size of matrix representing pupil}
  \item{obstruct}{central obstruction \emph{fraction}}
  \item{iname}{short string for identification}
}
\details{
  It's important to note that \code{zcoef} is treated differently than in other functions
  that use the same variable name. The first element \emph{must} be a piston (constant) term,
  which is stripped off and passed to \code{\link{pupil}} as the \code{piston} argument. The length
  of \code{zcoef} therefore should be one more than the length of \code{zlist$n}.
}
\value{
  A \code{size} by \code{size} matrix of intensity levels in the simulated interferogram.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  The relationship between wavefront phase and intensity is
  \code{iwf <- cos(2 * pi * wf + pi)}, which is the value returned by
  \code{synth.interferogram}. The plot routine in this function
  plots the image on a 256 level grayscale.
}
\seealso{
  \code{\link{Zernike}},
  \code{\link{pupil}}.
}
\examples{
zcoef <- c(0, 3, 3, rnorm(length(zlist.qf$n)-2, mean=0, sd = .01))

temp <- synth.interferogram(zcoef, iname="Random wavefront")

# lets see what it looks like in a star test

zcoef <- zcoef[-1]
zcoef[1:2] <- 0
fraunhofer(zcoef, displaywf=TRUE)
}
\keyword{Graphics}
\keyword{hplot}
