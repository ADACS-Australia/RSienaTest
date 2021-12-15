\name{bayesTest}
\alias{bayesTest}
\alias{testBayes}
\alias{simpleBayesTest}
\alias{multipleBayesTest}
\alias{print.multipleBayesTest}
\alias{plot.multipleBayesTest}
\title{
  Tests for sienaBayes results with print and plot methods
}
\description{
  These functions compute tests based on \code{\link{sienaBayesFit}} objects
  resulting from \code{\link{sienaBayes}}.
  Print and plot methods are available for the results of the
  multi-parameter test.
}
\usage{
simpleBayesTest(z, nfirst=z$nwarm+1, tested0=0,
                probs = c(0.025,0.975), ndigits=4)
multipleBayesTest(z, testedPar, nfirst=z$nwarm+1, tested0=0, ndigits=4)
\method{print}{multipleBayesTest}(x, descriptives=FALSE, ...)
\method{plot}{multipleBayesTest}(x, xlim=NULL, ylim=NULL, main=NULL, ...)
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{z}{A \code{\link{sienaBayesFit}} object,
    resulting from a call to \code{\link{sienaBayes}.}}
  \item{nfirst}{The first element of the MCMC chain used for calculating
      properties of the chain; this can be the first element for which
      it is assumed that convergence has occurred.}
  \item{tested0}{The value to be tested; for \code{simpleBayesTest} this
     must be a number (applied to all coordinates),
     for \code{multipleBayesTest} it can be a number
     (applied to all coordinates) or a vector.}
  \item{probs}{A vector of two numbers between 0 and 1, the credibility limits
  for the credibility intervals.}
  \item{ndigits}{Number of digits to be given when representing the limits
   of the credibility intervals.}
  \item{testedPar}{A vector with the numbers of the parameters
   to be tested (numbered as in the output of \code{simpleBayesTest}), or
   a \code{k * p} matrix, for testing \code{k} linear combinations
   given by the rows, where \code{p} is the number of effects in the model
   except the basic rate parameters.}
  \item{x}{A \code{\link{multipleBayesTest}} object,
    resulting from a call to \code{multipleBayesTest}.}
  \item{\dots}{Extra arguments for \code{print.multipleBayesTest}
     and \code{plot.multipleBayesTest}.}
  \item{descriptives}{Boolean: whether to print the posterior mean, s.d., and
     covariance matrix of the parameters included in the test.}
  \item{xlim}{Parameter \code{xlim} in \code{\link{plot}}; if \code{NULL},
   a sensible default will be used.}
  \item{ylim}{Parameter \code{ylim} in \code{\link{plot}}; if \code{NULL},
   a sensible default will be used.}
  \item{main}{Parameter \code{main} in \code{\link{plot}}; if \code{NULL},
    'posterior distances'.}
}
\details{
 Note that for the default values of \code{probs} and \code{tested0},
 \code{simpleBayesTest} gives information also contained in
 \code{\link{print.sienaBayesFit}}.\cr
 To use \code{multipleBayesTest}, it is advisable first to execute
 \code{simpleBayesTest} to see the numbers that index the \code{testedPar}
 parameter;  or its columns, if it is a matrix.\cr
 \code{multipleBayesTest} tests the hypothesis that the multivariate parameter
 for the effects defined by \code{testedPar} is equal to \code{tested0}.
 If  \code{testedPar} is a set or vector of numbers then the parameters
 with these numbers are tested; if \code{testedPar} is a \code{k * p} matrix,
 the tested hypothesis is for the \code{k} linear combinations defined by its
 rows, and the number of columns \code{p} should be equal to the number of
 non-rate effects,  which is the same as the number of effects listed by
 \code{simpleBayesTest}.\cr
 For the test, distances (for one parameter) or Mahalanobis distances
 (for more than one parameter) of the elements
 of the posterior sample from the posterior mean are calculated,
 and the \code{p}-value is the relative frequency that these are greater than the
 distance between the tested value and the posterior mean.\cr
 The \code{p}-values reported by \code{simpleBayesTest} are one-sided, and those
 reported by \code{multipleBayesTest} are two-sided for \code{k=1},
 and 'all-sided' if \code{k} is larger than 1.\cr
 The plot presents a density plot of the posterior Mahalanobis distances,
 and the observed distance is indicated by a vertical line (if within
 the plot window).
}
\value{\code{simpleBayesTest} produces a data frame containing for the non-rate
 effects their names,
 an indication of whether or not they are varying parameters,
 and test results.
 For each parameter separately, (for the varying parameters) the population
 mean or (for the non-varying parameters) the value is tested.
 The test results are the posterior credibility intervals for the
 probability limits defined by \code{probs},
 and the posterior probabilities of values larger than \code{tested0}.\cr
 \code{multipleBayesTest} produces an object of type \code{multipleBayesTest}
 which is a list with values:
 \item{prob}{posterior \code{p}-value}
 \item{chisquared}{test statistic}
  \item{postDistances}{sample from posterior of distances from \code{tested0}}
  \item{nullValue}{input value of \code{tested0}}
  \item{effectNames}{names of effects tested as given in \dots}
  \item{posteriorSample}{sample from posterior of tested effects}
 \code{print.multipleBayesTest} prints the posterior \code{p}-value with an
     explanation; and the proportions of all sign patterns
     (of deviations from 0) in the sample from the posterior.\cr
 \code{plot.multipleBayesTest} plots the density function of
    \code{x$postDistances}.\cr
  }
\references{See the manual and \url{http://www.stats.ox.ac.uk/~snijders/siena/}
}
\author{Tom Snijders}
\seealso{\code{\link{sienaBayes}}, \code{\link{print.sienaBayesFit}}}
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
  # Note this also shows the "randomEffects" column.
  FourAlgo <- sienaAlgorithmCreate(projname = "FourGroups", maxlike=TRUE)
\dontrun{
  bayes.model <- sienaBayes(FourAlgo, data = FourGroups,
        effects = FourEffects, nwarm=10, nmain=25, nrunMHBatches=10)
  bayes.model
  simpleBayesTest(bayes.model, probs = c(0.05,0.95))
  (mbt <- multipleBayesTest(bayes.model, c(2,3)))
  (mbt2 <- multipleBayesTest(bayes.model, c(2,3), tested0 = c(2,0.6)))
  plot(mbt)
  plot(mbt2)
  plot(mbt$posteriorSample)
  # Possible because the number of tested parameters is 2.
  A <- matrix(c(1,1,0),1,3)
  multipleBayesTest(bayes.model, A)
}
}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{ tests }