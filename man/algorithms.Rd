\name{algorithms}
\alias{algorithms}
\title{
  Function to estimate \code{Siena} model using alternative algorithms.
}
\description{
 Estimates a longitudinal network model using different methods from those
 in \code{\link{siena07}}.
}
\usage{
algorithms(data, effects, x, ...)
}
\arguments{
  \item{data}{
A \code{\link{siena}} data object}
  \item{effects}{
A \code{\link{sienaEffects}} object
}
  \item{x}{
A \code{\link{sienaModel}} object}

  \item{\dots}{
Other arguments controlling the algorithm
}
}
\value{
  A \code{\link{sienaFit}} object.
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}

\author{
Ruth Ripley}
\note{Not for general use}


\examples{
\dontrun{
mynet1 <- sienaDependent(array(c(tmp3, tmp4), dim=c(32,32,2)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
ans <- algorithms(mydata, myeff, sienaModelCreate())
}}
\keyword{ models }
