\name{interferogram}
\alias{interferogram}
\title{Interferogram object}
\description{
  Creates an instance of an interferogram object.
}
\usage{
interferogram(filename)
}
\arguments{
  \item{filename}{The name of a graphics file in \dQuote{Portable Anymap} format}
}
\details{
  Provides an \dQuote{object oriented} framework for analysis of interferograms.
  A call of the form \code{.thisint <- interferogram("myinterferogram.ppm")}
  creates an \dQuote{instance} of an interferogram object. All functions for processing
  the interferometry data are returned with the function call.

  The current version only supports graphics files in \dQuote{Portable Anymap} format,
  and RGB files are converted to grayscale using the R channel only on the assumption
  that almost all interferograms will use a He-Ne source.

  If \code{interferogram()} is called from the \code{\link{Rfringe}} GUI 
  \bold{ImageMagick's} (\url{http://www.imagemagick.org})
  \code{convert} tool is used if necessary to convert from any graphics format supported by \bold{ImageMagick} to
  pnm format. This means that \code{convert} must be present and on the search path if this feature
  is to be utilized. In Windows it also means that the path to the \bold{ImageMagick} installation must be placed
  earlier in the path list than system files, since there is also a system \code{convert} command.
}
\value{
    A list with the following components:
  \item{ev}{The environment of the \code{interferogram()} function call}
  \item{isInterferogram}{Identifies this object as an interferogram}
  \item{image.info}{Data entry function - Basic image information}
  \item{analysis.info}{Data entry function - Information required in wavefront analysis}
  \item{target.conic.info}{Data entry function - Target conic in single pass tests}
  \item{circle.pars}{Image analysis function - Outline aperture edge}
  \item{obstruct.pars}{Image analysis function - Outline obstruction (perforation)}
  \item{plot.fringes}{Plot function - Plots fringe trace}
  \item{autotrace}{Image analysis function - (Semi) automatic fringe tracing}
  \item{clearpoints}{Image analysis function - Fringe editing}
  \item{addpoints}{Image analysis function - Fringe editing}
  \item{clearfringe}{Image analysis function - Fringe editing}
  \item{retrace}{Image analysis function - Fringe editing}
  \item{addsegment}{Image analysis function - Fringe editing}
  \item{insertfringe}{Image analysis function - Fringe editing}
  \item{fitzernikes}{Image analysis function - Least squares fit to fringe centers}
  \item{plot.si}{Plot function - Synthetic interferogram}
  \item{plot.wf}{Plot function - Wavefront map}
  \item{plot.contour}{Plot function - Cross sections of wavefront along selected diameters}
  \item{plot.wf3d}{Plot function - Interactive 3d plot of wavefront (requires package \code{rgl})}
  \item{plot.residuals}{Plot function - Some possibly useful diagnostic plots of residuals from least squares fit}
  \item{plot.startest}{Plot function - Star test simulation with optional MTF plot}
  \item{print.summary}{Output function - Prints basic summary statistics to the console window}
  \item{print.details}{Output function - Prints more details of estimated Zernike coefficients}
  \item{print.latex}{Output function - A printable detailed report, in pdf format (requires \code{pdflatex})}
}
\note{
  These functions will rarely be accessed directly at the command line if the GUI is used.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}. Thanks to Steven Koehler for valuable programming
  ideas, especially his approach to object instantiation in \R.}
\seealso{
  \code{\link{Rfringe}},
  \code{\link{project}},
  \code{pixmap}.
}
\examples{
  \dontrun{
    .thisint <- interferogram("myinterferogram.ppm")  # create an interferogram object
    .thisint$circle.pars()	#interactively outline the aperture edge
    .thisint$obstruct.pars()	#outline the obstruction
    .thisint$autotrace()	#trace the fringes
    .thisint$fitzernikes()	#do the least squares fit
    .thisint$plot.wf()		#plot the wavefront
    .thisint$print.summary()	#print some basic summary stats
  }
}
\keyword{Statistics}
\keyword{models}
