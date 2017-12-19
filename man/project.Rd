\name{project}
\alias{project}
\title{Interfometry project}
\description{
  Creates an instance of an interferometry project.
}
\usage{
project(project.id, project.notes=NULL, project.tester=NULL, project.date=NULL)
}
\arguments{
  \item{project.id}{Character string - an identifier for the project}
  \item{project.notes}{Optional character string - additional project notes}
  \item{project.tester}{Optional character string - the tester (or any other useful information)}
  \item{project.date}{Optional character string - date of the test (or other information)}
}
\details{
  Provides an "object oriented" framework for the analysis of multiple
  interferograms as a group.
  A call of the form \code{.thisproject <- project("myproject")}
  creates an "instance" of a project object. All functions for processing
  the grouped data are returned with the function call.
}
\value{
    A list with the following components:
  \item{ev}{The environment of the \code{project()} function call}
  \item{isIntProject}{Identifies this object as an interferometry project}
  \item{project.addto}{Utility function - copies data from an interferogram into project object}
  \item{project.removefrom}{Utility function - removes interferogram data from project object}
  \item{plot.wf}{Plot function - Wavefront map}
  \item{plot.contour}{Plot function - Cross sections of wavefront along selected diameters}
  \item{plot.startest}{Plot function - Star test simulation with optional MTF plot}
  \item{plot.wf3d}{Plot function - Interactive 3d plot of wavefront (requires package \code{rgl})}
  \item{plot.spm}{Plot function - Coplots of RMS, P-V, and Strehl from individual interferograms}
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
  \code{\link{interferogram}},
}
\examples{
  \dontrun{
    .thisproject <- project("myproject")  # create a project object
    .thisproject$project.addto(.thisint$ev)	#note you add an interferogram 
    						# by passing its environment
    .thisproject$plot.wf()	#Plot averaged wavefront
    .thisproject$plot.spm()	#Possibly interesting comparisons of interferograms
    .thisproject$print.summary()	#Some basic summary statistics
  }
}
\keyword{Statistics}
\keyword{models}
