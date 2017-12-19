\name{fraunhofer}
\alias{fraunhofer}
\title{Star test simulator}
\description{
  A simple star test simulator, with optional MTF or wavefront plot.
}
\usage{
fraunhofer(zcoef, zlist = zlist.qf, obstruct = 0, lambda = 1, defocus = 5, 
    pupilsize = 255, npad=1024, gamma=2, psfmag = 2, 
    displaymtf = TRUE, displaywf = FALSE, fileout = FALSE)
}
\arguments{
  \item{zcoef}{Vector of Zernike coefficients}
  \item{zlist}{A list with named components \code{n, m, t} describing the contents of \code{zcoef}}
  \item{obstruct}{central obstruction \emph{fraction}}
  \item{lambda}{Wavelength. Defaults to 1 for a wavefront in waves.}
  \item{defocus}{Amount of defocus for side panes, in waves P-V.}
  \item{pupilsize}{Size of matrix representing pupil}
  \item{npad}{Size of the 0 padded matrix for fft calculations}
  \item{gamma}{gamma value for image display}
  \item{psfmag}{Magnification factor for in focus PSF display}
  \item{displaymtf}{Boolean: display MTF?}
  \item{displaywf}{Boolean: display Wavefront?}
  \item{fileout}{Boolean: Is the intended output a file?}
}
\details{
  Computes the intensity distribution from a monochromatic point source in the
  approximation of Fraunhofer diffraction theory. See \cite{Suiter (1994)} for
  a full discussion of the application of the theory to practical star testing
  of small telescopes.
}
\references{
  Born, M. and Wolf, E. 1999, \emph{Principles of Optics, 7th Edition},
  Cambridge University Press.

  Suiter, H. R., 1994, \emph{Star Testing Astronomical Telescopes},
  Willman-Bell, Inc.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  The ratio of npad to pupilsize determines the size of the Airy disk in the
  \dQuote{in focus} image. The default choices will produce an Airy disk with
  a diameter of approximately 10 pixels, and the amount of zero padding
  is sufficient to prevent significant aliasing for defocus values up to about
  20 waves.
  
  If memory limitations cause errors satisfactory images can still be obtained with,
  say, \code{npad=512} and \code{pupilsize=127}.
}
\seealso{
  \code{\link{Zernike}},
  \code{\link{pupil}}.
}
\examples{
# a random, but probably almost diffraction limited, wavefront

fraunhofer(rnorm(length(zlist.qf$n), mean=0, sd=0.01), displaywf=TRUE)
}
\keyword{Graphics}
\keyword{hplot}
