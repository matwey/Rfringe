\name{zlist.qf}
\alias{zlist.qf}
\title{Fringe set of Zernike Polynomials}
\description{List of the \sQuote{Fringe} set of Zernike polynomials}
\usage{
zlist.qf
}
\details{
  The \dQuote{Fringe} aka \dQuote{QuickFringe} set of Zernike polynomials.
}
\value{
  A list with the following components:
  \item{n}{Radial polynomial order}
  \item{m}{Azimuthal order}
  \item{t}{character for trig function: one of c("n", "c", "s")}
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  This is the default list of Zernike polynomial indices used
  by all functions that work with lists of Zernikes. It includes all
  Zernikes from 2nd through 10th order, plus the 12th order spherical
  term.
}
\keyword{Mathematics}
