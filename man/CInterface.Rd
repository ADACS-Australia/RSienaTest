\name{sienacpp}
\alias{sienacpp}
\keyword{models}
\title{Estimate parameters in a Siena model}
\usage{
sienacpp(x,
    nThreads=1,
    logLevelConsole='WARNING', logLevelFile='INFO',
    logBaseName=x$projname, logIncludeLocation=F,
    tt=NA, batch=NA, verbose=NA, silent=NA, initC=NA, useCluster=NA,
    clusterType=NA, clusterString=NA, nbrNodes=NA, clusterIter=NA,
    \dots)
}
\description{
Estimates parameters in a Siena model using Method of Moments,
or Generalized Method of Moments, based on direct simulation;
or using Maximum Likelihood by MCMC simulation.
Estimation is done using a Robbins-Monro algorithm.
Note that the data and particular model to be used must be passed
in using named arguments as \code{\dots}, and the specification for the
algorithm must be passed on as x, which is a \code{\link{sienaAlgorithm}}
object as produced by \code{\link{sienaAlgorithmCreate}}; see examples below.
}
\arguments{
\item{x}{
A control object, of class \code{\link{sienaAlgorithm}}.
}
\item{nThreads}{
Number of OpenMP threads each process can use.
}
\item{logLevelConsole}{
Minimal level of logging output to the console. One of 'DEBUG', 'VERBOSE',
'INFO', 'WARNING', 'ERROR', 'FATAL'. Where 'DEBUG' prints all messages,
'FATAL' only the messages with highest priority.
}
\item{logLevelFile}{
Like \code{logLevelConsole} but for logging to files.
}
\item{logBaseName}{
Base name of the logging files. Each process logs into a separate file. The
final log file name is \emph{<logBaseName>-<mpiRank>:<threadID>.log}.
}
\item{logIncludeLocation}{
When \code{TRUE} also log the code location (file, method, line) where the
logging message originated.
}
\item{\dots}{
Arguments for \code{\link{initializeFRAN}}. Usually, \code{data} and
\code{effects}, as in the examples below.
}
\item{batch}{Deprecated. See \code{logLevel}}
\item{verbose}{Deprecated. See \code{logLevel}}
\item{silent}{Deprecated. See \code{logLevel}}
\item{clusterString}{
Deprecated. The cluster setup is detected automatically via MPI.
}
\item{useCluster}{Deprecated. See \code{clusterString}}
\item{nbrNodes}{Deprecated. See \code{clusterString}}
\item{clusterIter}{Deprecated. See \code{clusterString}}
\item{clusterType}{Deprecated. See \code{clusterString}}
\item{initC}{Deprecated. C++ anyway.}
\item{tt}{Deprecated. No GUI supported.}
}
\details{
This implementation is an alternative to \code{\link{siena07}}.
In \code{sienacpp} the algorithm is coded entirely in C++,
whereas for \code{\link{siena07}} the 'front end' is coded in R.
The Generalized Method of Moments is available only in \code{sienacpp}.
Objects produced by \code{\link{siena07}} contain some more material
(usually not needed). For the rest, the two functions do the same.\cr

\code{sienacpp} runs a Robbins-Monro algorithm for parameter estimation using
for the Method of Moments the three-phase
implementation in Snijders (2001) and Snijders, Steglich and Schweinberger
(2007), with (if \code{findiff=TRUE} in the \code{\link{sienaAlgorithm}}
object) derivative estimation as in Schweinberger and Snijders (2007). Phase 1
does a few iterations to estimate the derivative matrix of the targets with
respect to the parameter vector. Phase 2 does the estimation. Phase 3 runs a
simulation to estimate standard errors and check convergence of the model.

The Generalized Method of Moments is used if any of the effects
specified in the effects object is defined with \code{type=gmm}.
For the Generalized Method of Moments the algorithm
of Amati, Schoenenberger, and Snijders (2015) is used,
which uses a three-phase Robbins-Monro algorithm as above,
but with an estimated optimal linear combination of the statistics.

If  \code{maxlike=TRUE} in the \code{\link{sienaAlgorithm}} object,
estimation is done by Maximum Likelihood implemented as in
Snijders, Koskinen and Schweinberger (2010) also using the three-phase
Robbins-Monro algorithm.

The deprecated parameters are included to allow using the same syntax
as for \code{\link{siena07}}.\cr
The \code{Dolby} option in \code{\link{sienaAlgorithmCreate}}
is not yet supported by \code{sienacpp}.\cr
\code{sienaFit} objects created by \code{sienacpp} cannot be used
for \code{\link{sienaTimeTest}} or \code{\link{sienaGOF}}.

The specification of the effects object for GMoM estimation requires
that in the effects object, apart from the basic rate effects,
some of the effects were specified
in  \code{\link{includeEffects}} with \code{type='eval'}
(the default, which means that this does not need to be stated)
and the others with \code{type='gmm'}.
The first then are evaluation effects defining the model specification,
the second are statistics used for estimation.
The method requires that the number of statistics (\code{type='gmm'})
is equal to or larger than the number of
evaluation effects (\code{type='eval'}).
See the example below.
}
\value{
Returns an object of class \code{\link{sienaFit}}, some parts of which are:
\item{OK}{Boolean indicating successful termination}
\item{termination}{Character string, values: "OK", "Error", or "UserInterrupt".
  "UserInterrupt" indicates that the user asked for early termination before
  phase 3.}
\item{f}{Various characteristics of the data and model definition.}
\item{theta}{Fitted value of theta.}
\item{pp}{Length of theta: number of parameters,
   excluding conditioning rate parameters (if any).}
\item{qq}{Number of statistics used in estimation: equal to \code{pp} for MoM
   and ML estimation.}
\item{covtheta}{Estimated covariance matrix of theta; this is not available if
  the sienaAlgorithm object x was produced with simOnly=TRUE.}
\item{dfra}{Matrix of estimated derivatives of expected statistics
(linear combinations used in Robbins-Monro updates) by parameters.}
\item{gamma}{Matrix of estimated derivatives of expected statistics (all)
by parameters.}
\item{sf}{Matrix of deviations from targets in phase 3.}
\item{sf2}{Array of statistics from simulations in phase 3.}
\item{tstat}{t-statistics for convergence.}
\item{targets}{Observed statistics.}
\item{targets2}{Observed statistics by wave, starting with the second wave .}
\item{ssc}{Score function contributions for each wave for each simulation in
  phase 3. Zero if finite difference method is used}
\item{Phase3nits}{Number of iterations actually performed in phase 3.}
Writes text output to the file named "projname.out", where projname is defined
in the \code{\link{sienaAlgorithm}} object \code{x}.
}
\author{Felix Schoenenberger, some modifications by Tom Snijders}
\seealso{\code{\link{sienaAlgorithmCreate}}, \code{\link{siena07}}.\cr
  There are print, summary and xtable methods for \code{\link{sienaFit}}
  objects: \code{\link{xtable}}, \code{\link{print.sienaFit}}.}
\references{
\itemize{
\item Amati, Viviana, Schoenenberger, Felix, and Snijders, Tom A.B. (2015).
Estimation of stochastic actor-oriented models for the evolution of
networks by generalized method of moments.
\emph{Journal de la Societe Francaise de Statistique}, 156, 140-165.
\item Schweinberger, Michael, and Snijders, Tom A.B. (2007). Markov models
for digraph panel data: Monte Carlo-based derivative estimation.
\emph{Computational Statistics and Data Analysis} 51, 4465-4483.
\item Snijders, Tom A.B. (2001). The statistical evaluation of social
network dynamics. \emph{Sociological Methodology}, 31, 361-395.
\item Snijders, Tom A.B., Steglich, Christian E.G., and Schweinberger,
Michael (2007). Modeling the co-evolution of networks and behavior.
Pp. 41-71 in \emph{Longitudinal models in the behavioral and related sciences},
edited by Kees van Montfort, Han Oud and Albert Satorra; Lawrence Erlbaum.
\item Steglich, Christian E. G., Snijders, Tom A. B., and Pearson, Michael
A. (2010). Dynamic networks and behavior: Separating selection from
influence. \emph{Sociological Methodology}, 40, 329-393.
\item Further see \url{http://www.stats.ox.ac.uk/~snijders/siena/ }
}}
\examples{
# Method of Moments:
algo <- sienaAlgorithmCreate(nsub=2, n3=100, seed=123)
net <- sienaDependent(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
dataset <- sienaDataCreate(net)
eff <- getEffects(dataset)
eff
(ans <- sienacpp(algo, data=dataset, effects=eff))
# or if a previous 'on track' result ans was obtained
(ans <- sienacpp(algo, data=dataset, effects=eff, prevAns=ans))
# Generalized Method of Moments:
algo <- sienaAlgorithmCreate(nsub=2, n3=100, dolby=FALSE, seed=123)
eff <- includeEffects(eff, density)
eff <- includeEffects(eff, density, type='gmm')
eff <- includeEffects(eff, recip)
eff <- includeEffects(eff, recip, realrecip, persistrecip, type='gmm')
eff
(ans <- sienacpp(algo, data=dataset, effects=eff))
}
\keyword{models}
