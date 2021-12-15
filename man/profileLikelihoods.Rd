\name{profileLikelihoods}
\alias{profileLikelihoods}
\title{
Function to calculate profile likelihood for a Siena model}
\description{
Calculates the profile likelihood for a Siena model for one or two parameters}
\usage{
    profileLikelihoods(resp, x, data, effects, i, j = NULL,
    gridl = c(0.8, 1.2), seqlen = 5, maxit = 2, method = "BFGS",
    trace = 0, nIter = 100, ...)
}
\arguments{
  \item{resp}{
	Value returned from a call to \code{\link{algorithms}} or a list
	  containing a set of parameter values, in a vector called \code{theta}.
}
  \item{x}{
A \code{\link{sienaAlgorithm}} object}
  \item{data}{
A \code{\link{siena}} data object
}
  \item{effects}{
A \code{\link{sienaEffects}} object
}
  \item{i}{
Required. Index of parameter to maximise over.}
  \item{j}{
Optional. Index of second parameter to maximise over.
}
\item{gridl}{
  Multiplier used to create grid for evaluation
}
\item{seqlen}{
  Number of points for evaluation.
}
\item{maxit}{
  Number of iterations of the optimizer \code{\link{optim}}.
}
  \item{method}{
Method for \code{\link{optim}}}
  \item{trace}{
trace argument for \code{\link{optim}}
}
  \item{nIter}{
Number of samples to use}
  \item{\dots}{
Other arguments to be passed to the initialize call.
}
}
\value{
  A \code{\link{sienaFit}} object, with added component \code{zz}
  containing the points used for the plot.
  Produces a plot.
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{
Ruth Ripley}

\note{Not for general use}

\seealso{
\code{\link{algorithms}}
}
\examples{
\dontrun{
mynet1 <- sienaDependent(array(c(tmp3, tmp4), dim=c(32,32,2)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
profileLikelihoods(list(theta=c(5,-2, 2)), sienaModelCreate(maxlike=TRUE),
   mydata, myeff, i=1)
}
}
\keyword{ models }
