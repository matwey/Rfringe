\name{Rfringe}
\alias{Rfringe}
\title{Rfringe GUI}
\description{
  GUI front end for R fringe analysis.
}
\usage{
Rfringe()
}
\arguments{
  none
}
\details{
  A Tcl/Tk wrapper for the functions provided in \code{\link{interferogram}}
  and \code{\link{project}}.
}
\note{
  The full functionality of \code{Rfringe} will be described in a separate user guide. This note briefly
  describes how objects are created and manipulated within the GUI interface.
  
  On startup \code{Rfringe} creates \code{NULL} valued objects named \code{.thisint} and \code{.thisproject}
  in the user's global environment. When the user subsequently creates a new interferogram object with the
  menu selection \code{File|Load interferogram from image...} \code{\link{make.names}} is run on the
  value entered in the image id field to create a new (syntactically correct) variable 
  in the user's workspace. This variable is
  assigned the value of the call to \code{interferogram()} and \code{.thisint} is in turn
  evaluated to the newly created interferogram object. Subsequent commands operate directly
  on \code{.thisint}.
  
  Similarly when a new project object is created with the menu selection \code{File|Project|Create...}
  a new variable is created from the Project ID entry, it is assigned the value of the call to
  \code{project()}, and \code{.thisproject} is evaluated to the newly created project object.
  
  The objects created with \code{make.names} are deliberately \emph{not} assigned unique names. This allows
  the user to replace an object simply by recycling her image or project id entries. The Rfringe GUI will
  warn if this is attempted.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}.
  GUI based in part on the package \code{Rcmdr} by John Fox.
}
\seealso{
  \code{tcltk},
  \code{\link{interferogram}},
  \code{\link{project}}
}
\examples{
# start the GUI

Rfringe()
}
\keyword{Statistics}
\keyword{models}
