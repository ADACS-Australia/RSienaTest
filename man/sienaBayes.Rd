\name{sienaBayes}
\alias{sienaBayes}
\alias{sienaBayesFit}
\alias{bayes}
\alias{glueBayes}
\concept{estimation}
\title{A function for fitting Bayesian models}
\description{A function to fit hierarchical Bayesian models random effects
  to sienaGroup data objects. Uses the function \code{\link{maxlikec}} for the
  SAOM part, the Bayesian part is performed in R.
}
\usage{
sienaBayes(data, effects, algo,  saveFreq=100,
    initgainGlobal=0.1, initgainGroupwise = 0.02, initfgain=0.2, gamma=0.05,
    initML=FALSE, priorSigEta=NULL,
    priorMu=NULL, priorSigma=NULL, priorDf=NULL, priorKappa=NULL,
    priorRatesFromData=2,
    frequentist=FALSE, incidentalBasicRates=FALSE,
    reductionFactor=0.5, delta=1e-04,
    nprewarm=50, nwarm=50, nmain=250, nrunMHBatches=20,
    nSampVarying=1, nSampConst=1, nSampRates=0,
    nImproveMH=100, targetMHProb=0.25,
    lengthPhase1=round(nmain/5), lengthPhase3=round(nmain/5),
    storeAll=FALSE, prevAns=NULL, usePrevOnly=TRUE,
    prevBayes=NULL, newProposalFromPrev=(prevBayes$nwarm >= 1),
    silentstart = TRUE,
    nbrNodes=1, clusterType=c("PSOCK", "FORK"),
    getDocumentation=FALSE)

glueBayes(z1, z2, nwarm2=0)
}
\arguments{
  \item{data}{A sienaGroup object as returned by \code{\link{sienaGroupCreate}}.
    It is planned to also allow a \code{\link{siena}} data object as returned by
    \code{\link{sienaDataCreate}.}}
  \item{effects}{sienaEffects object as returned by \code{\link{getEffects}(data)}.
    The effects indicated by its column \code{randomEffects} get
    effects that are varying across the groups.}
  \item{algo}{Algorithm object, as created by
    \code{\link{sienaAlgorithmCreate}}. Should contain all options required
    for the MCMC scheme, and a random seed if required.}
  \item{saveFreq}{Integer. If this is larger than 1, the provisional results
    are saved after each multiple of saveFreq iterations in the main phase,
    in a file with name PartialBayesResult.RData (if a file with this name
    exists, it will be overwritten).
    This file contains an object \code{z} of class \code{\link{sienaBayesFit}},
    with the provisional results.
    This is to guard against crashes or power failure.
    It can be used as value for \code{prevBayes} as indicated below.}
  \item{initgainGlobal}{Step sizes in initial searches for good parameter values
    across the groups.}
  \item{initgainGroupwise}{Step sizes in initial searches for good parameter
    values by group; can be up to 0.1 for larger networks,
    0 for very small networks.}
  \item{initfgain}{Positive number, used only for frequentist estimation
    and incidentalBasicRates:
    the gain factor in the Robbins Monro algorithm is \code{initfgain} during
    warming and Phase 1, and
    \code{initfgain * ((iteration number after Phase 1)^(-gamma))} after that.}
  \item{gamma}{Positive number, used only for frequentist estimation:
    see \code{initfgain}.}
  \item{initML}{Boolean, whether to use maximum likelihood estimation for the
     groupwise initial estimation.}
  \item{priorSigEta}{Vector of length equal to the number of fixed parameters,
     or \code{NULL}. If not \code{NULL}, values that are not \code{NA} are prior
     variances of \code{eta} (fixed effects); these must all be positive.}
  \item{priorMu}{Vector of length equal to the number of randomly varying
     parameters, or \code{NULL}.
     Prior mean of \code{mu} (global population mean for varying parameters);
     default: 0.
     For the basic rate parameters the defaults are data-dependent,
     but if ((\code{priorRatesFromData = 0}) and \code{priorMu=NULL}),
     for these parameters the prior mean will be 2.}
  \item{priorSigma}{Square matrix of dimension equal to length of \code{priorMu},
     or \code{NULL}.
     Prior global population covariance matrix \code{Sigma} for the varying parameters;
     default: identity matrix.}
  \item{priorDf}{Prior degrees of freedom for \code{Sigma} (global population
    covariance matrix); default: number of randomly varying parameters + 2.}
  \item{priorKappa}{Proportionality constant between prior covariance matrix
    and covariance matrix of prior distribution for \code{mu}; default: 1.}
  \item{priorRatesFromData}{-1, 0, 1, or 2. Determines the prior distribution
     for the rate parameters.\cr
     -1: all basic rate parameters are fixed; this must correspond to
         \code{effects};\cr
     0: prior is taken from \code{priorMu} and \code{priorSigma};\cr
     1: prior is defined by mean and covariance matrix of estimated rate parameters
        from initialization phase;\cr
     2: prior is defined by robust estimates for location and multivariate scale
        of estimated rate parameters from initialization phase.}
  \item{frequentist}{Currently only \code{frequentist=FALSE} works.
    Boolean: chooses between frequentist or Bayesian
    estimation of the global parameters.
    Frequentist estimation is possible only for at least 2 groups.}
  \item{incidentalBasicRates}{Boolean. If this is \code{TRUE}, the basic rate
    parameters are defined specifically for each group, and estimated
    using a Robbins Monro algorithm; if \code{FALSE}, they have a
    common prior for all groups.}
  \item{reductionFactor}{Positive number.
   If \code{priorRatesFromData} = 1 or 2, the prior distribution for the rate
   parameters is the estimated covariance matrix of the estimated rate parameters
   in the initialization phase multiplied by \code{reductionFactor},
   plus a contribution of 0.01 for the variances.
   For many data sets the default will be OK. If the rate parameters diverge,
   as evidenced by the traceplots of the posterior values,
   then a smaller value is advisable, such as 0.01.
   Has no effect if \code{priorRatesFromData} <= 0.}
  \item{delta}{When the global population covariance matrix becomes
    non-positive definite (i.e., has one or more negative correlations)
    during iterations, it is changed so that
   all correlations are at least delta.}
  \item{nprewarm}{Number of iterations in the pre-warm-up phase,
    before \code{improveMH}.}
  \item{nwarm}{Number of iterations in the warm-up phase.
    Used only if \code{prevBayes} is \code{NULL}.
    Then it should be at least 5.}
  \item{nmain}{Number of iterations in the main phase.
     Should be at least 10.}
  \item{nrunMHBatches}{Integer: thinning ratio in MCMC process;
     but thinning at the lowest level is further determined by parameter
     \code{mult} in \code{algo},
     see \code{\link{sienaAlgorithmCreate}}. Should be even.}
  \item{nSampVarying}{Number of samples of varying parameters
        for each chain sample.}
  \item{nSampConst}{Number of samples of constant parameters ("\code{eta}")
        for each chain sample.}
  \item{nSampRates}{Number of extra samples of basic rate parameters
        for each chain sample.}
  \item{nImproveMH}{Number of iterations per \code{improveMH} step.
        If \code{nImproveMH=0}, no \code{improveMH} steps at all.}
  \item{targetMHProb}{Desired proportion of acceptances in MH steps.
       This can be one number, or a vector of two numbers. In the latter case,
       the first number applies to the MH steps per group (for each group),
       the second to the MH steps for the constant parameters ("\code{eta}").}
  \item{lengthPhase1}{Only used for frequentist estimation or
      \code{incidentalBasicRates}:
      length of the first phase of the Robbins Monro algorithm.
      \code{lengthPhase1 + lengthPhase3} should be strictly 
	  less than \code{nmain}.}
   \item{lengthPhase3}{Only used for frequentist estimation or
      \code{incidentalBasicRates}:
      length of the third phase of the Robbins Monro algorithm.
      \code{lengthPhase1 + lengthPhase3} should be strictly 
	  less than \code{nmain}.}
  \item{storeAll}{Boolean: whether to store parameters for all MCMC iterations,
      i.e., before thinning. \code{storeAll=TRUE} may lead to producing
      very large objects and is not recommended for usual operation.}
  \item{prevAns}{An object of class "sienaFit" as returned by
      \code{\link{siena07}} for the same data set and effects object.
      This then is the result of a multi-group estimation for these data,
      from which scaling information (derivative matrix and standard deviation of
      the deviations) can be extracted along with the parameters estimates
      which will be used as the initial values, unless
      \code{algo} requests the use of standard initial values.\cr
      If \code{prevAns=NULL}, then a multi-group estimation
      using \code{\link{siena07}} will be performed as part of the
      initialization. If \code{prevAns} is an object as described,
      estimated with unconditional estimation, with 500 or more runs in Phase 3,
      then if \code{usePrevOnly=TRUE} this multi-group estimation will be skipped;
      if \code{usePrevOnly=FALSE} a multi-group estimation still will be
      performed, but starting at the results provided by \code{prevAns};
      if \code{usePrevOnly=TRUE} but with less than 500 runs in Phase 3,
      then the multi-group estimation will be done only for Phase 3
      (skipping Phases 1 and 2). }
  \item{usePrevOnly}{Boolean: see immediately above.}
  \item{prevBayes}{An object of class \code{sienaBayes} as returned by function
      \code{sienaBayes}, on which the current function will continue.
      For example, the object \code{z} contained in the intermediate saved result
      \code{PartialBayesResult.RData} as mentioned above.
      Initialization and warming phases are skipped. \cr
      If these are given, the values of \code{nrunMHBatches},
      \code{nSampVarying}, \code{nSampConst}, and \code{nSampRates} supersede
      those in \code{prevBayes}.
      If these are given, the values of \code{prevAns}, \code{nwarm},
      \code{incidentalBasicRates}, \code{priorRatesFromData}, and the prior
      distributions are disregarded.\cr
      A new set of \code{nmain} iterations in the main phase are done.
      All further relevant parameters of the function call should be identical.
      This is not completely checked, so errors
      may occur if this is not the case.}
  \item{newProposalFromPrev}{Boolean, has consequences only when
      \code{prevBayes} is given: whether to use a robust estimate of the
      covariance matrix of parameters in \code{prevBayes} to define a new
      proposal distribution.  Else the proposal distribution is the same
      as in the \code{prevBayes} object.}
  \item{silentstart}{Boolean: whether to suppress most information
      to the console during the calculation of initial values.}
  \item{nbrNodes}{Number of processes to be used. Cannot be more than
      the number of periods summed over number of groups.}
   \item{clusterType}{If using multiple processes, whether to use
      forking processes or not. (Only "PSOCK" can be used on Windows.)}
   \item{getDocumentation}{Flag to allow documentation of internal
      functions, not for use by users.}
   \item{z1}{\code{sienaBayes} object.}
   \item{z2}{\code{sienaBayes} object with the same data and model specification
      as \code{z1}.}
   \item{nwarm2}{Number of warming iterations in \code{z2}; the first
      \code{nwarm2} iterations of \code{z2} will be left out of the
      combined object.}
}
\details{
  The function \code{sienaBayes} is for Bayesian estimation of one group or
  of multiple groups all having the same number of waves and the same
  model specification.
  It wraps Bayesian sampling of parameters around calls to
  \code{\link{maxlikec}}.  The RSiena manual has a lot of information. \cr
  Effects can be randomly varying between groups, or global (i.e., constant
  across groups).  This is indicated by the keyword \code{random} in the
  call of \code{\link{setEffect}}, and reported when using the
  keyword \code{includeRandoms} in \code{\link{print.sienaEffects}}.\cr
  For the groupwise parameters normal distributions are assumed with conjugate
  priors.
  The conjugate prior is an inverse Wishart distribution for the covariance
  matrix \code{Sigma}, with parameters \code{Lambda^{-1}} and \code{priorDf},
  where \code{Lambda = priorDf*priorSigma};
  and, conditional on \code{Sigma}, for the expected value a multivariate normal with
  mean \code{priorMu} and covariance matrix \code{Sigma/priorKappa}.\cr
  For the fixed parameters, if \code{priorSigEta=NULL}, the improper constant
  prior is used.
  Else, for the components of \code{eta} (fixed effects parameter) for which
  \code{priorEta} is \code{NA},
  the improper constant prior is used, and for other components,
  independent normal priors with mean 0 and variances given by \code{priorEta}.
  \cr
  Recommendation: use the normal prior only for effects of group-level
  variables.  If a non-zero prior mean is desired, transform the
  group-level variable.  \cr
  The required dimensions of the prior parameters \code{priorSigEta},
  \code{priorMu} and \code{priorSigma} depend on the number of groupwise varying
  parameters and are given in the output for \code{\link{print.sienaEffects}}
  when using the keyword \code{includeRandoms} .\cr
  If \code{priorRatesFromData} is 1 or 2,
  the prior distribution for the basic rate parameters is determined
  in a data-dependent way.  For the non-varying parameters,
  a flat prior is assumed. \cr
  The frequentist option currently is not supported (probably broken).\cr
  The procedure consists of three parts: initialization, warming,
  main phase.\cr
  In the initialization phase, initial parameter values and the proposal
  covariance matrix for Metropolis-Hastings steps for groupwise parameters
  are obtained from, first, Method of Moments estimation of a parameter vector
  assumed to be the same across the groups (in a multi-group estimation);
  replacing this estimate by a precision-weighted mean of this estimate
  and the prior mean;
  followed by one subphase of the Robbins-Monro algorithm for Method of
  Moments estimates for the groups separately, starting from the overall
  estimate, with step size \code{initgain}.
  This is skipped if \code{initgain=0}.\cr
  If the \code{prevBayes} object is supplied, this initialization phase is
  skipped, and if \code{newProposalFromPrev} the proposal covariance matrices
  are calculated from the simulated chains in this object.\cr
  From this basis, \code{nprewarm} Metropolis Hastings steps are taken;
  however, if no acceptances occur in two consecutive steps, the prewarming
  phase is ended.
  The proposal covariance matrices then are scaled, in the
  function '\code{improveMH}', to achieve a proportion of accepted MH proposals
  approximately equal to \code{target}.\cr
  After initialization and scaling of the proposal covariance matrices,
  a warming phase is done of \code{nwarm} Bayesian proposals
  each with a number of MH steps, followed again by the function '\code{improveMH}'.
  If the \code{prevBayes} object is supplied, the warming phase is skipped.
  \cr
  Finally \code{nmain} repeats of \code{nrunMHBatches} are performed of a
  number of MH steps sampling chains, plus \code{nSampVarying} MH steps
  sampling the varying parameters ('\code{theta_j}') plus \code{nSampConst}
  MH steps sampling the non-varying parameters ('\code{eta}') plus one Gibbs
  step sampling the global mean and covariance matrix of the varying parameters
  ('\code{mu}' and '\code{Sigma}').
  In the warming as well as the final phase, the number of MH steps
  within each run is determined by parameter \code{mult}
  ("multiplication factor") in the algorithm object \code{algo}.

  The function is time-consuming.  When starting to use it, it is advisable
  to start with low values of \code{nmain} to explore computing time.\cr
  Initialisation is an important part, and will be shorter if \code{prevAns}
  is given.  If an earlier \code{sienaBayes} object is available for
  this data specification (which may differ in which effects are specified
  as being random), then if the name of this object is \code{zz}
  it will be possible to use \code{prevAns=zz$initialResults}.

  When the procedure seems to diverge, and for very small groups, it is
  advisable to use smaller values of the parameters \code{initgainGlobal}
  and \code{initgainGroupwise}.

  \code{glueBayes(z1,z2)} combines two existing \code{sienaBayesFit} objects
  with the same data and model specifications (this is not entirely
  checked) into one by putting \code{z2} after \code{z1}.
  This is intended to be used when \code{z2} was produced
  by calling \code{sienaBayes} with \code{prevBayes=z1}.

  For the function \code{\link{print.sienaEffects}}, the options
  \code{includeRandoms=TRUE} and \code{dropRates=TRUE} are meant to be useful
  especially for effects objects used for \code{sienaBayes}.
}
\value{
  \code{sienaBayes} as well as \code{glueBayes} return an object of class
  \code{sienaBayesFit}. This is a list containing, among other things:
  \item{priorMu}{prior global population mean (not quite the same
  as corresponding input parameter)}
  \item{priorSigma}{prior global population covariance matrix
  (not quite the same as corresponding input parameter)}
  \item{priorKappa}{proportionality constant between prior covariance matrix
                and covariance matrix of prior distribution for the mean}
  \item{priorDf}{prior degrees of freedom for covariance matrix}
  \item{nmain}{Number of iterations used in the main phase, as given
     in the call of \code{sienaBayes}}
  \item{nwarm}{Number of iterations used in the warm-up phase, as given
     in the call of \code{sienaBayes}; 0 if \code{prevBayes} was used.
     The value of \code{nwarm} is used in other
     functions such as \code{\link{print.sienaBayesFit}},
     and when the object is used for \code{sienaBayes} as a
     \code{prevBayes} parameter with \code{newProposalFromPrev=TRUE}.
     If it was deemed that converge begins later (or earlier) than
     the value of \code{nwarm}, this value can be changed for further
     use of the object produced by \code{sienaBayes}.}
  \item{effectName}{array of names of effects included in the model}
  \item{f$groupNames}{array of names of groups included in the model}
  \item{initialResults}{sienaFit object: result of initial MoM estimation
  under the assumption of same parameters across groups}
  \item{ThinParameters}{array of dimensions (\code{nwarm+nmain} iterations
    by parameters by groups): sampled groupwise parameters}
  \item{ThinPosteriorMu}{array of dimensions (\code{nwarm+nmain} iterations
    by parameters): sampled global means of varying parameters}
  \item{ThinPosteriorEta}{array of dimensions (\code{nwarm+nmain} iterations
    by fixed parameters): sampled fixed parameters}
  \item{ThinPosteriorSigma}{array of dimensions (\code{nwarm+nmain} iterations
    by parameters by parameters): sampled global covariance matrix
    of varying parameters}
  \item{nrunMH}{vector (by groups * waves) of number of MH steps in the
   innermost loop of the likelihood simulations of the SAOM;
   this is a data-dependent multiple of the value of \code{mult}
   set in \code{\link{sienaAlgorithmCreate}}}
  \item{ThinBayesAcceptances}{matrix of number of acceptances
   in the \code{nrunMHBatches} MH steps in each iteration,
   and summed over dependent variables,
   for each group and (duplicated) for the vector of fixed effects}
  \item{acceptances}{if \code{storeAll=TRUE}: matrix of booleans:
   whether the corresponding change to the parameters was accepted
   in each step, by group}
  \item{MHacceptances}{if \code{storeAll=TRUE}:
   array of acceptances of the MH steps, by step type and group
   but summed over dependent variables}
  \item{MHrejections}{if \code{storeAll=TRUE}:
   array of rejections of the MH steps}
  \item{MHproportions}{if \code{storeAll=TRUE}:
   array of proportions of the MH steps accepted}
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}

  Koskinen, J.H. and T.A.B. Snijders (2007).
  Bayesian inference for dynamic social network data.
  \emph{Journal of Statistical Planning and Inference}, 13, 3930-3938.
}
\author{Ruth Ripley, Johan Koskinen, Tom Snijders }
\seealso{\code{\link{siena07}}, \code{\link{sienaGroupCreate}},
   \code{\link{bayesTest}}, \code{\link{extract.sienaBayes}},
   \code{\link{plotPostMeansMDS}}\cr
  There are print and summary methods for sienaBayesFit objects, see
  \code{\link{print.sienaBayesFit}}.}
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
  print(bayes.model, nfirst=20)
  summary(bayes.model)
  bayes.nextModel <- sienaBayes(FourAlgo, data = FourGroups,
       effects = FourEffects, nmain=15, nrunMHBatches=10,
       prevBayes = bayes.model)
  bayes.combinedModel <- glueBayes(bayes.model, bayes.nextModel)
  summary(bayes.combinedModel)
}
}
\keyword{models}
