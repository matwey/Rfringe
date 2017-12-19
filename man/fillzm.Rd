\name{fillzm}
\alias{fillzm}
\title{Create a matrix of Zernike polynomial values}
\description{
  Creates a matrix of Zernike polynomials values from vectors \code{rho}
  and \code{theta}
}
\usage{
fillzm(rho, theta, phi = 0, zlist = zlist.qf)
}
\arguments{
  \item{rho}{normalized radius, \eqn{0 <= rho <= 1}}
  \item{theta}{angular coordinate}
  \item{phi}{angular coordinate to rotate entire coordinate system}
  \item{zlist}{A list with named components \code{n, m, t}}
}
\details{
  The arguments \code{rho} and \code{theta} must be vectors of the same length.

  The optional argument \code{phi} is an angle specified in \emph{degrees}
  for rotation of the entire coordinate system. All angles are measured
  increasing \emph{counterclockwise}.
}
\value{
  A matrix of size \code{length(rho)} rows by \code{length(zlist$n)} columns.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  The main use for this function in Rfringe is creating a matrix of Zernikes to
  be used as the predictors in the least squares fit to fringe orders.
}
\seealso{
  \code{\link{Zernike}},
  \code{\link{makezlist}},
  \code{\link{zlist.qf}}.
}
\keyword{Mathematics}
