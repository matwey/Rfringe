\name{listInterferograms}
\alias{listInterferograms}
\alias{listProjects}
\title{List Interferogram and Project Objects}
\description{
  Lists the objects identified as Interferogram or Project
  objects in the user's workspace.
}
\usage{
listInterferograms()
listProjects()
}
\value{
  Character vectors containing the names of objects identified as
  interferograms or interferogram projects.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  Interferogram objects are identified by having a named attribute \code{isInterferogram}.
  
  Similarly, project objects have a named attribute \code{isIntProject}.
}
\seealso{
  \code{\link{interferogram}},
  \code{\link{project}}.
}
\keyword{Basics}
\keyword{list}

