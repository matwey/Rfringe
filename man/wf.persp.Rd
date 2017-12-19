\name{wf.persp}
\alias{wf.persp}
\title{3D Wavefront plot}
\description{
  An alternative 3D Wavefront plot using the \R base package plotting function
  \code{\link{persp}}.
}
\usage{
wf.persp(wf, zoom.wf = 1, theta=0, phi=30, ...)
}
\arguments{
  \item{wf}{A wavefront matrix as returned by \code{\link{pupil}}}
  \item{zoom.wf}{Zoom factor to stretch wavefront heights}
  \item{theta}{Value of theta to pass to \code{\link{persp}}}
  \item{phi}{Value of phi to pass to \code{\link{persp}}}
  \item{\dots}{Additional parameters for call to \code{persp}}
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\seealso{
  \code{\link{pupil}},
  \code{\link{persp}},
  \code{\link{wf.3dplot}}.
}
\example{
# yet another random wavefront

wf <- pupil(zcoef=rnorm(length(zlist.qf$n), mean=0, sd=.01))

wf.persp(wf)
}
\keyword{Graphics}
\keyword{hplot}
