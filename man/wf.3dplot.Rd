\name{wf.3dplot}
\alias{wf.3dplot}
\title{Interactive 3D Wavefront plot}
\description{
  Uses the \code{\link[pkg:rgl]{rgl}} package for OpenGL graphics to
  produce an interactive 3D wavefront map.
}
\usage{
wf.3dplot(wf, zoom.wf = 1)
}
\arguments{
  \item{wf}{A wavefront matrix as returned by \code{\link{pupil}}}
  \item{zoom.wf}{Zoom factor to stretch wavefront heights}
}
\details{
  The \code{rgl} package is available at \url{http://wsopuppenkiste.wiso.uni-goettingen.de/~dadler/rgl/}.
  The README file in the source distribution contains installation instructions.
}
\references{
\url{http://wsopuppenkiste.wiso.uni-goettingen.de/~dadler/rgl/}
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\seealso{
  \code{\link{pupil}}
}
\example{
# make sure we've got the library
require(rgl)

# yet another random wavefront

wf <- pupil(zcoef=rnorm(length(zlist.qf$n), mean=0, sd=.01))

wf.3dplot(wf)
image(wf, col=topo.colors(256), asp=1)
contour(wf, add=TRUE)
}
\keyword{Graphics}
\keyword{interactive}
