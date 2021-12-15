\name{extract.sienaBayes}
\alias{extract.sienaBayes}
\alias{extract.Bayes}
\alias{extract.posteriorMeans}
\title{
  Extraction of posterior samples or posterior means from sienaBayes results
}
\description{
  The first function extracts posterior samples from a list of
  \code{\link{sienaBayesFit}} object to be used, e.g.,
  for assessing convergence.\cr
  The second function extracts posterior means and standard deviations per group
  from a \code{\link{sienaBayesFit}} object.
}
\usage{
extract.sienaBayes(zlist, nfirst=zlist[[1]]$nwarm+1, extracted,
                   sdLog=TRUE)
extract.posteriorMeans(z, nfirst=z$nwarm+1, pmonly=1, excludeRates=FALSE,
                   verbose=TRUE)
}
\arguments{
  \item{zlist}{A list of \code{\link{sienaBayesFit}} objects,
    further called 'chains',
    resulting from calls to \code{\link{sienaBayes}} with a common
    data set and model specification.}
  \item{z}{A \code{\link{sienaBayesFit}} object.}
  \item{nfirst}{Integer: the first element for the first MCMC chain
      used for calculating properties of the chain.}
  \item{extracted}{The parameters for which posterior samples are to be
     extracted: \cr
     \code{"all"}: all parameters;\cr
     \code{"rates"}: all groupwise rate parameters;\cr
     \code{"varying"}: all varying non-rate parameters:
                      global means and standard deviations;\cr
     \code{"non-varying"}: all estimated non-varying
                   (and therefore, non-rate) parameters;\cr
     \code{"objective"}: all non-rate parameters.
     }
  \item{sdLog}{Logical: should the logarithms be taken of the posterior
      standard deviations (applies only to varying parameters).}
  \item{pmonly}{if \code{pmonly=0} the posterior means and standard deviations
       are extracted, if \code{pmonly=1} only the posterior means,
      if \code{pmonly=2} only the posterior standard deviations.}
  \item{excludeRates}{Logical: whether to exclude the rate parameters.}
  \item{verbose}{Logical: gives some console output to show ongoing activity
       (may be reassuring for large number of groups).}
} 
\details{
 \code{extract.sienaBayes} produces a 3-dimensional array, iterations by chains
 by parameters, suitable for use, e.g., in function \code{monitor()}
 of package \code{rstan}. The 'iterations' are draws from the
 posterior distribution (provided that convergence was achieved),
 after the thinning in \code{sienaBayes} implied by its parameters
 \code{nrunMHBatches}, \code{nSampVarying}, \code{nSampConst} and
 \code{nSampRate}.\cr
 The number of iterations should be the same for all chains.\cr
 It is not checked that the specifications of the chains are identical.
 For the parameter names, the shortNames in the first chain are used.

 \code{extract.posteriorMeans} produces the posterior means of the
 groupwise varying parameters. The resulting matrix has the groups
 in the rows. The effects are in the columns: for each effect comes first
 the posterior mean, then the posterior standard deviation.
}
\value{For \code{extract.sienaBayes} a 3-dimensional array, iterations
  by chains by parameters,
  the third dimension having the names of the parameters.\cr
  For \code{extract.posteriorMeans} a matrix with the groups in the rows
  and all effects in the columns, with for each effect two columns:
  first the posterior mean (effect name preceded by "p.m.") and then the
  posterior standard deviation (effect name preceded by "psd.").\cr
  If some names are duplicated (e.g., because of the presence of
  evaluation and creation effects), these can be changed by the user
  after creation of the array.
  }
\references{See the manual and \url{http://www.stats.ox.ac.uk/~snijders/siena/}
}
\author{Tom Snijders}
\seealso{\code{\link{sienaBayes}, \link{plotPostMeansMDS}}}
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
  bayes.models <- replicate(4, sienaBayes(FourAlgo, data = FourGroups,
        effects = FourEffects, nwarm=10, nmain=25, nrunMHBatches=10),
        simplify=FALSE)
  bayes.extracted <- extract.sienaBayes(bayes.models, extracted="all")
  dim(bayes.extracted)
  dimnames(bayes.extracted)
  extract.posteriorMeans(bayes.models[[1]])
}
}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{ tests }