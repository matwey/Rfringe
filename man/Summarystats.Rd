\name{Summarystats}
\alias{pupilrms}
\alias{pupilpv}
\alias{strehlratio}
\title{Wavefront summaries}
\description{
  Estimate the RMS or P-V wavefront error over the pupil given in the argument.
}
\usage{
pupilrms(pupil)
pupilpv(pupil)
strehlratio(rms)
}
\arguments{
  \item{pupil}{\code{pupil} is the matrix created by the call to \code{\link{pupil}}}
  \item{rms}{The rms wavefront error}
}
\value{
  An estimate of the RMS or P-V error of the wavefront, or Strehl ratio.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  The function \code{pupilrms} simply returns the standard deviation of the defined values in
  \code{pupil}, which is a crude but usually good enough approximation to
  the properly defined integral over the aperture.

  For a wavefront defined entirely in terms of a vector \code{zcoef} of Zernike coefficients
  \code{rms <- sqrt(crossprod(zcoef))} is faster and more accurate.

  \code{pupilpv} does the obvious. There is no analytical solution in general for P-V.

  \code{strehlratio} computes Mahajan's approximation to the Strehl ratio.
}
\seealso{
  \code{\link{pupil}}.
}
\examples{
# A random vector of Zernike coefficients

zcoef <- rnorm(length(zlist.qf$n), mean=0, sd=0.01)
wf <- pupil(zcoef=zcoef)
image(wf, col=topo.colors(256), asp=1)
contour(wf, add=TRUE)
pupilrms(wf)
sqrt(crossprod(zcoef)) #should be the same to about 4 digits
pupilpv(wf)
strehlratio(sqrt(crossprod(zcoef))) #probably around 0.8
}
\keyword{Mathematics}
