\name{print.sienaBayesFit}
\alias{print.sienaBayesFit}
\alias{summary.sienaBayesFit}
\alias{print.summary.sienaBayesFit}
\alias{shortBayesResults}
\alias{plotPostMeansMDS}
\title{Methods and functions for processing sienaBayes objects}
\usage{
\method{print}{sienaBayesFit}(x, nfirst=NULL, ...)

\method{summary}{sienaBayesFit}(object, nfirst=NULL, ...)

\method{print}{summary.sienaBayesFit}(x, nfirst=NULL, ...)

shortBayesResults(x, nfirst=NULL)

plotPostMeansMDS(x, pmonly=1, excludeRates=TRUE, nfirst=NULL, ...)
}
\description{
  \code{print} and \code{summary} methods for \code{\link{sienaBayesFit}} objects,
  and further functions for interpretation of results.
}
\arguments{
  \item{object}{An object of class \code{\link{sienaBayesFit}}.}
  \item{x}{An object of class  \code{\link{sienaBayesFit}} or
    \code{summary.sienaBayesFit} as appropriate.}
  \item{nfirst}{The first element of the MCMC chain used for calculating
      properties of the chain; i.e., the first element for which
      it is assumed that convergence has occurred. If \code{NULL},
      the value \code{x$nwarm+1} will be used.}
  \item{pmonly}{if \code{pmonly=0} a plot is made for posterior means and standard deviations,
      if \code{pmonly=1} only for the posterior means,
      if \code{pmonly=2} only for the posterior standard deviations.}
  \item{excludeRates}{Logical: whether to exclude the rate parameters
      for determining the coordinates.}
  \item{\dots}{For extra arguments. At present used only for
     further specification of the plot made by \code{plotPostMeansMDS},
    allowing further arguments for \code{\link{plot}} and \code{\link{text}}.}
}
\value{
  The function \code{print.sienaBayesFit} prints results (estimated posterior means,
  posterior standard deviations, and percentiles) of a Bayesian analysis, as resulting
  from function \code{\link{sienaBayes}}.\cr
  The function \code{summary.sienaBayesFit} prints more extensive results,
  including groupwise results.\cr
  The columns "cred. from" and "cred. to" give the left and right end points
  of the estimated 95 percent credibility interval for the parameters, and the
  column "p" gives the estimated posterior probability that the parameter is
  greater than 0.\cr
  The function \code{shortBayesResults} produces a data frame containing
  posterior means and standard deviations, and 95 percent percentiles.
  This is meant to be a component for use in other functions.\cr
  Of this data frame, variables "name" to "effectName" are as in a
  \code{\link{sienaEffects}} object;
  variables "postMeanGlobal", "postSdGlobal", "cFrom", "cTo" are the posterior
  mean, posterior standard deviation, and lower and upper boundaries of the
  95 percent posterior interval for the global parameters (\code{mu} and
  \code{eta});
  variables  "postSdBetween", "cSdFrom", "cSdTo" are the posterior mean,
  and similar percentiles, for the global between-group standard deviations
  (\code{sigma}) of the parameters that are randomly varying between groups.

  The function \code{plotPostMeansMDS} makes an MDS plot of the posterior means.
  This can be useful for detection of outliers.
  Its return value is a list with two components:
  \code{corpm}, the similarities between the groups (correlation matrix of the
  transposed matrix of posterior means), used as input for the MDS;
  \code{points}, the matrix of two-dimensional MDS coordinates.

  It is advisable to employ the \code{first} parameter; looking at trace plots
  will be helpful to estimate the number of runs from which convergence
  has occurred; this may well be different from the default value \code{x$nwarm+1}.

  Function \code{\link{siena.table}} can also be used to make html or LaTeX
  tables of posterior means and standard deviations of a
  \code{\link{sienaBayesFit}} object.
}
\seealso{\code{\link{siena.table}}, \code{\link{sienaBayes}},
    \code{\link{extract.sienaBayes}}, \code{\link{extract.posteriorMeans}}.}
\references{
  See also \url{http://www.stats.ox.ac.uk/~snijders/siena/}
  \cr
  For the MDS method used in \code{plotPostMeansMDS}: see p. 308 of\cr
  W.N. Venables and B.D. Ripley (2002),
  \emph{Modern Applied Statistics with S.} Fourth edition. Springer.
}
\author{Tom Snijders}
\keyword{print}
\examples{
  Group1 <- sienaDependent(array(c(N3401, HN3401), dim=c(45, 45, 2)))
  Group3 <- sienaDependent(array(c(N3403, HN3403), dim=c(37, 37, 2)))
  Group4 <- sienaDependent(array(c(N3404, HN3404), dim=c(33, 33, 2)))
  Group6 <- sienaDependent(array(c(N3406, HN3406), dim=c(36, 36, 2)))
  dataset.1 <- sienaDataCreate(Friends = Group1)
  dataset.3 <- sienaDataCreate(Friends = Group3)
  dataset.4 <- sienaDataCreate(Friends = Group4)
  dataset.6 <- sienaDataCreate(Friends = Group6)
  FourGroups <- sienaGroupCreate(
        list(dataset.1, dataset.3, dataset.4, dataset.6))
  FourEffects <- getEffects(FourGroups)
  FourEffects <- includeEffects(FourEffects, transTrip)
  FourEffects <- setEffect(FourEffects, density, random=TRUE)
  FourEffects <- setEffect(FourEffects, recip, random=TRUE)
  print(FourEffects, includeRandoms=TRUE)
# Note this also shows the "randomEffects" column
# and the dimensions of the objects for specifying the priors.
  print(FourEffects, includeRandoms=TRUE, dropRates=TRUE)
# Note this does not show the rate effects.
  FourAlgo <- sienaAlgorithmCreate(projname = "FourGroups", maxlike=TRUE,
                                   seed=123)
\dontrun{
  bayes.model <- sienaBayes(FourAlgo, data = FourGroups,
       effects = FourEffects, nprewarm=10, nwarm=10, nmain=25, nrunMHBatches=10)
  bayes.model
# If one would assume that converge was reached by iteration 20:
  print(bayes.model, nfirst=20)
  summary(bayes.model)
  shortBayesResults(bayes.model)
# In the following, "main" is used as ... :
  plotPostMeansMDS(bayes.model, main="MDS of posterior means")
  siena.table(bayes.model, type="html", d=3)
}
}
