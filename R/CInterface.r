###############################################################################
# SIENA: Simulation Investigation for Empirical Network Analysis
#
# Web: http://www.stats.ox.ac.uk/~snijders/siena/
#
# File: CInterface.r
#
# Description: Provides an R interface to the C objects via external pointers.
###############################################################################

lsmem <- function(name=NULL, ...) {
  if (is.null(name)) {
    names  <- ls(...)
    napply <- function(fn) sapply(names, function(x) fn(get(x, ...)))
  } else {
    list <- name
    names(list) <- names(list)
    if (is.list(list)) {
      napply <- function(fn) sapply(list, function(x) fn(x))
    } else {
      return();
    }
  }
  obj       <- list()
  obj$class <- napply(function(x) as.character(class(x))[1])
  obj$mode  <- napply(mode)
  obj$size  <- napply(object.size)
  obj$hsize <- napply(function(x) capture.output(print(object.size(x), units="auto")))
  obj$len   <- napply(length)

  sum <- sum(obj$size)
  class(sum) <- "object_size"

  obj$class <- c(obj$class, "--")
  obj$mode  <- c(obj$mode, "--")
  obj$size  <- c(obj$size, sum)
  obj$hsize <- c(obj$hsize, capture.output(print(sum, units="auto")))
  obj$len   <- c(obj$len, sum(obj$len))

  data.frame(obj)[order(obj$size, decreasing=T),]
}

# Helper function to get the location (function name) of a dispatched log
# message.
# Mostly from: http://stackoverflow.com/questions/7307987/logging-current-function-name
clean.callstack <- function(callstack) {
  val <- sapply(callstack, function(call) {
                call <- strsplit(paste(call, collapse="\t"), "\t")[[1]]
                switch(call[1],
                       "lapply" = call[3],
                       "sapply" = call[3],
                       "do.call" = call[2],
                       "function" = "FUN",
                       "source" = "###",
                       "eval.with.vis" = "###",
                       call[1])
})
  val[grepl("\\<function\\>", val)] <- "FUN"
  #val <- val[!grepl("(###|FUN)", val)]
  val <- head(val, -1) # ignore the LOG function
  val <- tail(val, 3)  # print the last 3 in the stack
  paste(val, collapse=".")
}

# Wrapper for the C++ logging.
#   priority  Name of the priority level. Or 'R' to cat to console.
#   ...       Concatenated into the log message.
LOG <- function(priority, ..., callstack=sys.calls()){
  msg <- paste(sapply(list(...), function(x) {
                      if (is.character(x)) return(x);
                      if (length(x)==1) as.character(x);
                      paste(capture.output(print(x)), collapse="\n")
                       }), collapse=" ")
  if (priority == 'R') {
    #cat('[', clean.callstack(callstack), '] ', msg, '\n', sep='')
    cat(msg, '\n', sep='')
  } else {
    .Call(C_sienaLog, PACKAGE=pkgname,
			priority, paste(clean.callstack(callstack), msg, collapse="\n"))
  }
}

# Helper function to setup the C++ logger for the old siena07 functions.
sienaSetupLogger <- function(logLevelConsole='INFO', logLevelFile='DEBUG',
                             logBaseName, logIncludeLocation=F) {
  .Call(C_sienaSetupLogger, PACKAGE=pkgname,
		logLevelConsole, logLevelFile, logBaseName, logIncludeLocation, 1)
}

# The RNG setup use by siena07.
setupSeed <- function(x, z) {
  randomseed2 <-  NULL
  if (!is.null(x$randomSeed)) {
    set.seed(x$randomSeed, kind="default")
#    seed <- x$randomSeed
  } else {
    if (!nzchar(Sys.getenv("RSIENA_TESTING"))) {
      if (exists(".Random.seed")) {
        rm(.Random.seed, pos=1)
      }
      newseed <- trunc(runif(1) * 1000000)
      set.seed(newseed)  ## get R to create a random number seed for me.
#      seed <- NULL
    } else {
      newseed <- NULL
#      seed <- NULL
    }
  }
  z$randomseed2 <- randomseed2
}

# Sets up the RNG state equal to the parallel package setup.
# Mostly copied from the parallel package.
# See: R sources (R/src/library/parallel/R/RngStream.R)
mpiClusterSetRNGStream <- function(iseed) {
  if (.Call(C_sienaMPISize, PACKAGE=pkgname) > 1) {
    RNGkind("L'Ecuyer-CMRG")
    # Fix (as in the value can't change anymore) iseed (use it in some way,
    # print(iseed) would also work).
    set.seed(iseed)
    # Compare the following 3 snippets in case your wondering what 'fix' means
    # ##########                       # R 3.0.2     R 2.15.2
    #    f <- function(x) {
    #      print(x)                    # 570175512,  570175512
    #      RNGkind("L'Ecuyer-CMRG")
    #      print(x)                    # 570175512,  570175512
    #    }
    #    RNGkind("default")
    #    set.seed(1)
    #    f(as.integer(runif(1, max=.Machine$integer.max)))
    # ##########                       # R 3.0.2      R 2.15.2
    #    f <- function(x) {
    #      RNGkind("L'Ecuyer-CMRG")
    #      print(x)                    # 1658026146, 1658026146
    #    }
    #    RNGkind("default")
    #    set.seed(1)
    #    f(as.integer(runif(1, max=.Machine$integer.max)))
    # ##########                       # R 3.0.2     R 2.15.2
    #    f <- function(x) {
    #      require(parallel)
    #      RNGkind("L'Ecuyer-CMRG")
    #      print(x)                    # 1262921772, 1658026146
    #    }
    #    RNGkind("default")
    #    set.seed(1)
    #    f(as.integer(runif(1, max=.Machine$integer.max)))
    # ##########
    # require(parallel) # Not really, just to set up the seeds (nextRNGStream)
    set.seed(iseed)   # Now set the correct seed
    rank <- .Call(C_sienaMPIRank, PACKAGE=pkgname)
    seeds <- vector("list", rank+1)
    seeds[[1]] <- .Random.seed
    for(i in seq_len(rank)) seeds[[i+1]] <- parallel::nextRNGStream(seeds[[i]])
    LOG("VERBOSE", "Set cluster seed", iseed, seeds[[rank+1]])
    assign(".Random.seed", seeds[[rank+1]], envir=.GlobalEnv)
  } else {
    LOG("VERBOSE", "Use default seed.")
  }
}

# This is 'siena07' with the call to robmon replaced with the new one.
# See: man/CInterface.Rd for some documentation.
sienacpp <- function(x, # the thing returned from sienaModelCreate
                     nThreads=1,
                     # logger
                     logLevelConsole='WARNING', logLevelFile='INFO',
                     logBaseName=x$projname,
                     logIncludeLocation=F,
                     # deprecated parameters, still accepted but ignored
                     tt=NA, batch=NA, verbose=NA, silent=NA, initC=NA,
                     useCluster=NA, clusterType=NA, clusterString=NA,
                     nbrNodes=NA, clusterIter=NA,
                     # passed to initializeFRAN
                     ...
                     )
{
  rank <- .Call(C_sienaMPIRank, PACKAGE=pkgname)
  .Call(C_sienaSetupLogger, logLevelConsole, logLevelFile,
        paste(logBaseName, rank, sep='-'),
        logIncludeLocation, nThreads, PACKAGE=pkgname)

  # Log warnings when using deprecated parameters.
  if (!is.na(tt))
    LOG('WARNING', "'tt' is deprecated. No GUI supported.")
  if (!is.na(batch))
    LOG('WARNING', "'batch' is deprecated.",
        "Use 'logLevel' to control the amount of output.")
  if (!is.na(verbose))
    LOG('WARNING', "'verbose' is deprecated.",
        "Use 'logLevel' to control the amount of output.")
  if (!is.na(silent))
    LOG('WARNING', "'silent' is deprecated.",
        "Use 'logLevel' to control the amount of output.")
  if (!is.na(initC))
    LOG('WARNING', "'initC' is deprecated.",
        "C++ code should now initialize it self.")
  if (!is.na(useCluster) || !is.na(clusterType) || !is.na(clusterString)
      || !is.na(nbrNodes) || !is.na(clusterIter)) {
    LOG('WARNING', "'useCluster', 'clusterType', 'clusterString' and",
        "'nbrNodes' are deprecated. MPI setup is detected automatically.")
  }

  on.exit(function() {
          Report(closefiles=TRUE)
          RNGkind("default")
        })

  z <- NULL
  z$x <- x
  z$maxlike <- x$maxlike # initializeFRAN accesses both

  setupSeed(x, z)
  mpiClusterSetRNGStream(as.integer(runif(1, max=.Machine$integer.max)))

  is.batch(TRUE)
  z$pb <- list(pb=NULL, pbval=0, pbmax=1)
  Report(openfiles=TRUE, projname=x$projname)
  ## reset the globals for interrupts
  NullChecks()

  z <- initializeC(z, x, initC=FALSE, ...)
#browser()
  f <- FRANstore()
  z$gmm <- any(z$effects$type=='gmm')
  if (z$gmm) {
    ngmm <- length(which(z$effects$type == 'gmm'))
    neval <- length(which(z$effects$type == 'eval'))
    if (ngmm < neval) {
      LOG('WARNING', "Less 'gmm' type effects than 'eval' type effects. (", ngmm, "<", neval, ")")
      LOG('WARNING', "This will not work, forgot to add 'mandatory' effects to the gmm?")
      return(NULL)
    }
    if (x$dolby) {
      LOG('WARNING', "'gmm' does not work together with 'dolby'.")
      return(NULL)
    }
	if (any(z$effects$type=='creation') || (any(z$effects$type=='endow')))
	{
      LOG('WARNING', "'gmm' does not allow 'creation' or 'endow' effects.")
      return(NULL)
	 }
  }

  # Run the estimation.
  z$estimationtime <- proc.time()['elapsed']
  z$sienafit <- .Call(C_sienaEstimateGroup, PACKAGE=pkgname, z,
                      f$pModel, f$pData, nThreads)
  z$estimationtime <- proc.time()['elapsed'] - z$estimationtime
  z <- reformatSienaFit(z)
  # Since we have not set all fields terminateFRAN() would fail.
  terminateC(z)

  Report(closefiles=TRUE) # also mentioned in on.exit, but seems not to happen
  # LOG("DEBUG", "after terminateFRAN")
  # .Call("sienaPrintMemUsage", PACKAGE=pkgname)
  # LOG("DEBUG", "\ngc\n", paste(capture.output(gc()), collapse='\n'))
  # LOG("DEBUG", "\n.GlobalEnv\n", lsmem(pos=".GlobalEnv"))
  # LOG("DEBUG", "\nFRAN\n", lsmem(FRANstore()))
  # LOG("DEBUG", "\nx\n", lsmem(x))
  # LOG("DEBUG", "\nz\n", lsmem(z))

  if (rank == 0) {
    z$OK <- TRUE
    z$termination <- 'OK'
    class(z) <- "sienaFit"
    if (z$gmm) {
      LOG('INFO', "'sienaFit' object is from a 'gmm' run.")
    }
    return(z)
  } else {
    LOG('INFO', 'This is a slave process. Returning NULL.')
    return(NULL)
  }
}

# Copy and reformat the results form the first data set to match the format
# returned by RSiena. Doing this allows to run the existing printing methods.
reformatSienaFit <- function(z) {
  f <- z$sienafit[[1]]
  z$n1 <- f$n1[[length(z$sienafit[[1]]$n1)]]
  z$x$n3 <- f$n3
  z$n <- z$n1 + sum(f$n2) + z$n3 # total number of simulations
  z$theta <- f$theta
  z$tstat <- f$tstat  # superseded below for gmm
  z$rate <- f$rate
  z$vrate <- f$rateError
  z$fixed <- f$fixed
  z$gamma <- f$gamma
  z$gmmweight <- f$phase3_gmmweight
  z$dfra <- f$phase3_derivative
#  z$dfra <-  f$gamma
  z$covtheta <- f$phase3_covariance
  z$sf2 <- f$phase3_period_wise_statistics
  z$sf <- f$phase3_statistics
  z$ssc <- f$phase3_period_wise_scores
  z$covtheta <- f$phase3_covariance
  z$pp <- length(z$theta)
  z$qq <- dim(z$gamma)[1]
  z$msf <- as.matrix(cov(z$sf)) # as.matrix() just in case z$pp = 1
  z$se <- diag(as.matrix(z$covtheta))
  z$Phase3nits <- f$n3
  if (!z$withPrevans){
    z$targets <- f$targets
	z$targets2 <- t(z$sienafit[[1]]$period_wise_targets)
  }
#  z$sf <- t(apply(z$sf2, 1, function(x) colSums(x)) - z$targets) is the same

# TS: The following two lines show a way to get elements from
# the sienafit components; but since the multigroup option seems
# not implemented correctly, I leave the rest for now.
  z$rate <- unlist(lapply(z$sienafit, function(x){x$rate}))
  z$vrate <- unlist(lapply(z$sienafit, function(x){x$rateError}))

  if (length(z$sienafit) >= 2)
  {
	message('\nFor multigroup data sets, sienacpp probably is incorrect.')
  }

# TS: for computing tconv.max, in siena07
# there are further conditions, search for 'toosmall'; probably superfluous,
#  dmsf <- diag(z$msf)
#  sf <- colMeans(z$sf)
#  toosmall <- (dmsf < 1e-20 * z$scale * z$scale)
#  toosmall2 <- (abs(sf) < 1e-10 * z$scale)
#  tstat <- rep(NA, z$pp)
#  tstat[!toosmall]<- sf[!toosmall] / sqrt(dmsf[!toosmall])
#  tstat[toosmall & toosmall2] <- 0
#  tstat[toosmall & !toosmall2] <- NA
#  z$tstat <- tstat
  # tconv.max = Maximum value of t-ratio for convergence,
  # for any linear combination.
  z$tconv.max <- NA
  if (z$gmm)
  {
    W <- solve(z$msf * f$blockStructure)
	# if f$blockStructure is not identically 1, check that the following still is OK!
    B0 <- t(z$gamma) %*% W
    B1 <- solve(diag(sqrt(rowSums(B0*B0)))) %*% B0
	mean.dev <- colSums(z$sf %*% t(B1))/dim(z$sf)[1]
	cov.dev <- B1 %*% z$msf %*% t(B1)
	evalstat <- match(z$requestedEffects$functionName[z$requestedEffects$type
	                   %in% c('eval','rate')],
		z$requestedEffects$functionName[z$requestedEffects$type %in% c('gmm','rate')])
	z$tstat <- colMeans(z$sf[,evalstat])/sqrt(diag(z$msf[evalstat,evalstat]))
  }
  else
  {
	if (sum(!z$fixed) > 0)
	{
		mean.dev <- colSums(z$sf)[!z$fixed]/dim(z$sf)[1]
		cov.dev <- z$msf[!z$fixed,!z$fixed]
	}
	else
	{
		mean.dev <- 0
		cov.dev <- 1
	}
  }


	if (inherits(try(thisproduct <- solve(cov.dev, mean.dev), silent=TRUE),
				"try-error"))
			{
	Report('Overall maximum t-ratio for convergence not computable.\n', outf)
			}
	else
	{
		z$tconv.max <- sqrt(t(mean.dev) %*% thisproduct)
	}

  # Set (dim)names to effect names.
  # Simulation (parameter bound) effects (rate, eval, creation, endowment).
  paramNames <- z$requestedEffects$shortName[which(z$requestedEffects$type!='gmm')]
  # Statistics (gmm function) effects (rate, gmm).
  if (z$gmm) {
    statNames <- z$requestedEffects$shortName[
      which(z$requestedEffects$type=='gmm' | z$requestedEffects$type=='rate')]
  } else {
    statNames <- paramNames
  }
  if (length(z$theta) == length(paramNames)) {
    names(z$theta) <- paramNames
    names(z$fixed) <- paramNames
  }
  if (nrow(z$dfra) == length(statNames) && ncol(z$dfra) == length(paramNames)) {
    dimnames(z$dfra) <- list(statNames, paramNames)
  }
  if (nrow(z$gamma) == length(statNames) && ncol(z$gamma) == length(paramNames)) {
    dimnames(z$gamma) <- list(statNames, paramNames)
  }
  z
}

# Small version of terminateFRAN. Only clears the C++ pointers.
terminateC <- function(z) {
  f <- FRANstore()
  f$pModel <- NULL
  f$pData <- NULL
  z$f$myeffects <- lapply(z$f$myeffects, function(x){x$effectPtr <- NULL;x})
  FRANstore(NULL)
}
# Smaller version of initializeFRAN.
initializeC <- function(z, x, data, effects, prevAns=NULL, initC,
                        profileData=FALSE, returnDeps=FALSE,
                        returnChains=FALSE, byGroup=FALSE,
                        returnDataFrame=FALSE, byWave=FALSE,
                        returnLoglik=FALSE, onlyLoglik=FALSE)
{
  z$effectsName <- deparse(substitute(effects))
  ## fix up the interface so can call from outside robmon framework
  if (is.null(z$FinDiff.method)) z$FinDiff.method <- FALSE
  if (is.null(z$int)) z$int <- 1
  if (is.null(z$int2)) z$int2 <- 1
  # if (!initC) ## ie first time round
  # {
  if (!inherits(data,"siena")) stop("not valid siena data object")
  ## check the effects object
  defaultEffects <- getEffects(data)
  if (is.null(effects)) {
    effects <- defaultEffects
  } else {
    ## todo check that the effects match the data dependent variables
    userlist <- apply(effects[effects$include,], 1, function(x)
                      paste(x[c("name", "effectName",
                                "type", "groupName")],
                            collapse="|"))
    deflist <- apply(defaultEffects, 1, function(x)
                     paste(x[c("name", "effectName",
                               "type", "groupName")],
                           collapse="|"))
    if (!all(userlist %in% deflist)) {
      bad <- which(!(userlist %in% deflist))
      print(userlist[bad])
      stop("invalid effect requested: see above ")
    }
  }
  if (!inherits(effects, "data.frame")) {
    stop("effects is not a data.frame")
  }
  if (x$useStdInits) {
    if (any(effects$effectName != defaultEffects$effectName)) {
      stop("Cannot use standard initialisation with a ",
           "different effect list")
    }
    effects$initialValue <- defaultEffects$initialValue
  }
  ## get data object into group format to save coping with two
  ## different formats
  if (inherits(data, "sienaGroup")) {
    nGroup <- length(data)
  } else {
    nGroup <- 1
    data <- sienaGroupCreate(list(data), singleOK=TRUE)
  }
  ## add any effects needed for time dummies
  tmp <- sienaTimeFix(effects, data)
  data <- tmp$data
  effects <- tmp$effects
  if (!x$useStdInits) {
    if (!is.null(prevAns) && inherits(prevAns, "sienaFit")) {
      effects <- updateTheta(effects, prevAns)
    }
  }
  ## find any effects not included which are needed for interactions
  tmpEffects <- effects[effects$include, ]
  interactionNos <- unique(c(tmpEffects$effect1, tmpEffects$effect2,
                             tmpEffects$effect3))
  interactionNos <- interactionNos[interactionNos > 0]
  interactions <- effects$effectNumber %in%
  interactionNos
  effects$requested <- effects$include
  requestedEffects <- effects[effects$include, ]

  effects$include[interactions] <- TRUE
  effects <- effects[effects$include, ]

  ## split and rejoin both versions before continuing
  depvarnames <- names(data[[1]]$depvars)

  effects1 <- split(requestedEffects, requestedEffects$name)
  effects1order <- match(depvarnames, names(effects1))
  requestedEffects <- do.call(rbind, effects1[effects1order])
  row.names(requestedEffects) <- 1:nrow(requestedEffects)

  effects1 <- split(effects, effects$name)
  effects1order <- match(depvarnames, names(effects1))
  effects <- do.call(rbind, effects1[effects1order])
  row.names(effects) <- 1:nrow(effects)

  ## now set up z provisionally
  z$theta <- requestedEffects$initialValue
  z$fixed <- requestedEffects$fix
  z$test <- requestedEffects$test
  z$pp <- length(z$theta)
  z$posj <- rep(FALSE, z$pp)
  z$posj[requestedEffects$basicRate] <- TRUE
  z$BasicRateFunction <- z$posj
  z$scale <- rep(0.1, z$pp)

  ## sort out names of user specified interaction effects
  effects <- fixUpEffectNames(effects)
  ## copy interaction names to the requested effects
  requestedEffects$effectName <- effects[effects$requested,
                                         "effectName"]
  requestedEffects$functionName <- effects[effects$requested,
                                           "functionName"]

  if (attr(data, "compositionChange")) {
    if (x$maxlike) {
      stop("Not possible to use maximum likelihood estimation ",
           "with composition change")
    }
  }

  ## if not specified whether conditional or nor, set to conditional
  ## iff there is only one dependent variable (therefore number 1)
  ## and not maxlike
  if (is.na(x$cconditional)) {
    x$cconditional <- !x$maxlike && (length(depvarnames) == 1)
    if (x$cconditional) {
      x$condvarno <- 1
    }
  }
  types <- sapply(data[[1]]$depvars, function(x) attr(x, "type"))
  ## now check if conditional estimation is OK and copy to z if so
  z$cconditional <- FALSE
  if ((x$cconditional) & (!(attr(data, "compositionChange")))) {
    if (x$maxlike) {
      stop("Conditional estimation is not possible with",
           "maximum likelihood method")
    }
    ##  if (nets == 1) not sure if this is necessary
    ##  {
    z$cconditional <- TRUE
    ## find the conditioning variable
    observations <- attr(data, "observations")
    ## this is actual number of waves to process
    if (x$condname != "") {
      z$condvarno <- match(x$condname, attr(data, "netnames"))
      z$condname <- x$condname
    } else {
      z$condvarno <- x$condvarno
      z$condname <- attr(data, "netnames")[x$condvarno]
    }
    z$condtype <- attr(data, "types")[z$condvarno]
    if (z$condtype == "oneMode")
      z$symmetric  <-  attr(data, "symmetric")[[z$condvarno]]
    else
      z$symmetric <- FALSE
    ## find the positions of basic rate effects for this network
    z$condvar <- (1:nrow(requestedEffects))[requestedEffects$name==
                                 z$condname][1:observations]
    z$theta<- z$theta[-z$condvar]
    z$fixed<- z$fixed[-z$condvar]
    z$test<- z$test[-z$condvar]
    z$pp<- z$pp - length(z$condvar)
    z$scale<- z$scale[-z$condvar]
    z$BasicRateFunction <- z$posj[-z$condvar]
    z$posj <- z$posj[-z$condvar]
    z$theta[z$posj] <- z$theta[z$posj] / requestedEffects$initialValue[z$condvar]
  }

  ## unpack data and put onto f anything we may need next time round.
  f <- lapply(data, function(xx, x) unpackData(xx, x), x=x)

  attr(f, "netnames") <- attr(data, "netnames")
  attr(f, "symmetric") <- attr(data, "symmetric")
  attr(f, "allUpOnly") <- attr(data, "allUpOnly")
  attr(f, "allDownOnly") <- attr(data, "allDownOnly")
  attr(f, "allHigher") <- attr(data, "allHigher")
  attr(f, "allDisjoint") <- attr(data, "allDisjoint")
  attr(f, "allAtLeastOne") <- attr(data, "allAtLeastOne")
  attr(f, "anyUpOnly") <- attr(data, "anyUpOnly")
  attr(f, "anyDownOnly") <- attr(data, "anyDownOnly")
  attr(f, "anyHigher") <- attr(data, "anyHigher")
  attr(f, "anyDisjoint") <- attr(data, "anyDisjoint")
  attr(f, "anyAtLeastOne") <- attr(data, "anyAtLeastOne")
  attr(f, "types") <- attr(data, "types")
  attr(f, "observations") <- attr(data, "observations")
  attr(f, "compositionChange") <- attr(data, "compositionChange")
  attr(f, "exooptions") <- attr(data, "exooptions")
  attr(f, "groupPeriods") <- attr(data, "groupPeriods")
  attr(f, "periodNos") <- attr(data, "periodNos")
  attr(f, "numberNonMissingNetwork") <-
    attr(data, "numberNonMissingNetwork")
  attr(f, "numberMissingNetwork") <- attr(data, "numberMissingNetwork")
  attr(f, "numberNonMissingBehavior") <-
    attr(data, "numberNonMissingBehavior")
  attr(f, "numberMissingBehavior") <- attr(data, "numberMissingBehavior")

  ##  attr(f, "totalMissings") <- attr(data, "totalMissings")

  if (x$maxlike && x$FinDiff.method) {
    stop("Finite difference method for derivatives not available",
         "with Maximum likelihood method")
  }
  ## if any networks symmetric must use finite differences and not maxlike
  ## scores are now available (30.10.10) but still not maxlike.
  ## check model type: default for symmetric is type 2 (forcing model).
  syms <- attr(data,"symmetric")[ attr(data,"types") %in% c("oneMode","bipartite")]
  z$FinDiffBecauseSymmetric <- FALSE
  z$modelType <- x$modelType
#  if (any(!is.na(syms) & syms)) {
    ##     z$FinDiff.method <- TRUE
    ##     z$FinDiffBecauseSymmetric <- TRUE
    # if (x$maxlike) {
      #                stop("Maximum likelihood method not implemented",
      #                     "for symmetric networks")
    # }
  z$modelType[(z$modelType == 1) & syms] <- 2
  z$behModelType <- x$behModelType
#  }
  if (z$cconditional) {
    attr(f, "change") <- sapply(f, function(xx)attr(xx$depvars[[z$condname]],
                                 "distance"))
    attr(f,"condEffects") <- requestedEffects[z$condvar,]
    effcondvar <- (1:nrow(effects))[effects$name==
                        z$condname][1:observations]
    effects <- effects[-effcondvar, ]
    requestedEffects <- requestedEffects[-z$condvar,]
  }
  ## use previous dfra only if everything matches including method
  z$haveDfra <- FALSE
  z$withPrevans <- FALSE
  if (!is.null(prevAns) && inherits(prevAns, "sienaFit") &&
      prevAns$maxlike == z$maxlike && !is.null(prevAns$sienafit)) {
    nonGMMEfffects = requestedEffects[which(requestedEffects$type!='gmm'),]
    if (!is.null(prevAns$dfra) &&
        ncol(prevAns$gamma) == nrow(nonGMMEfffects) &&
        all(colnames(prevAns$gamma) == nonGMMEfffects$shortName) &&
        all(prevAns$fix == nonGMMEfffects$fix) &&
        !is.null(prevAns$sf)) {
      z$haveDfra <- TRUE
      z$gamma <- prevAns$gamma
      z$gmmweight <- prevAns$gmmweight
      z$dfra <- prevAns$dfra
      z$fixed <- prevAns$fixed
      z$sf <- prevAns$sf
	  z$targets <- prevAns$targets
	  z$targets2 <- prevAns$targets2
	  z$withPrevans <- TRUE
      z$dinv <- prevAns$dinv
      # z$dinvv must not be taken from prevAns,
      # because the value of diagonalize
      # is defined in x and may have changed.
      # Therefore here we copy the corresponding lines
      # from phase1.r.
      if (!x$diagg) {
        # Partial diagonalization of derivative matrix
        # for use if 0 < x$diagonalize < 1.
        temp <- (1-x$diagonalize)*z$dfra +
					x$diagonalize*diag(diag(z$dfra))
        temp[z$fixed, ] <- 0.0
        temp[, z$fixed] <- 0.0
        diag(temp)[z$fixed] <- 1.0
        # Invert this matrix
        z$dinvv <- solve(temp)
      }
      # check for backward compatibility with pre-1.1-220 versions:
      if (is.null(prevAns$regrCoef)) {
        z$regrCoef <- rep(0, z$pp)
        z$regrCor <- rep(0, z$pp)
      } else {
        z$regrCoef <- prevAns$regrCoef
        z$regrCor <- prevAns$regrCor
      }
    }
  }
  z$effects <- effects
  z$requestedEffects <- requestedEffects
  # }
  # else ## initC, i.e just send already set up data into new processes
  # {
  #   f <- FRANstore()
  #   ## Would like f to be just the data objects plus the attributes
  #   ## but need the effects later. Also a few other things,
  #   ## which probably could be attributes but are not!
  #   ## They will be automatically removed: if needed they must be readded
  #   ff <- f
  #   nGroup <- f$nGroup
  #   f[(nGroup + 1): length(f)] <- NULL
  # }
  pData <- .Call(C_setupData, PACKAGE=pkgname,
                 lapply(f, function(x)(as.integer(x$observations))),
                 lapply(f, function(x)(x$nodeSets)))
  ans <- .Call(C_OneMode, PACKAGE=pkgname,
               pData, lapply(f, function(x)x$nets))
  ans <- .Call(C_Bipartite, PACKAGE=pkgname,
               pData, lapply(f, function(x)x$bipartites))
  ans <- .Call(C_Behavior, PACKAGE=pkgname,
               pData, lapply(f, function(x)x$behavs))
  ans <-.Call(C_ConstantCovariates, PACKAGE=pkgname,
              pData, lapply(f, function(x)x$cCovars))
  ans <-.Call(C_ChangingCovariates, PACKAGE=pkgname,
              pData, lapply(f, function(x)x$vCovars))
  ans <-.Call(C_DyadicCovariates, PACKAGE=pkgname,
              pData, lapply(f, function(x)x$dycCovars))
  ans <-.Call(C_ChangingDyadicCovariates, PACKAGE=pkgname,
              pData, lapply(f, function(x)x$dyvCovars))
  ans <-.Call(C_ExogEvent, PACKAGE=pkgname,
              pData, lapply(f, function(x)x$exog))
  ## split the names of the constraints
  higher <- attr(f, "allHigher")
  disjoint <- attr(f, "allDisjoint")
  atLeastOne <- attr(f, "allAtLeastOne")
  froms <- sapply(strsplit(names(higher), ","), function(x) x[1])
  tos <- sapply(strsplit(names(higher), ","), function(x) x[2])
  ans <- .Call(C_Constraints, PACKAGE = pkgname, pData, froms[higher], tos[higher],
    froms[disjoint], tos[disjoint], froms[atLeastOne], tos[atLeastOne])

  ##store the address
  f$pData <- pData
  ## register a finalizer
  ans <- reg.finalizer(f$pData, clearData, onexit = FALSE)

  # if (!initC) {
    storage.mode(effects$parm) <- "integer"
    storage.mode(effects$group) <- "integer"
    storage.mode(effects$period) <- "integer"
    effects$effectPtr <- rep(NA, nrow(effects))
    splitFactor <- factor(effects$name, levels=attr(f, "netnames"))
    if (!all(attr(f,"netnames") %in% effects$name)) {
      stop("Must have at least one effect for each dependent variable")
    }
    myeffects <- split(effects, splitFactor)
    myCompleteEffects <- myeffects
    ## remove interaction effects and save till later
    basicEffects <-
      lapply(myeffects, function(x) {
             x[!x$shortName %in% c("unspInt", "behUnspInt"), ]
               }
    )
    basicEffectsl <-
      lapply(myeffects, function(x) {
             !x$shortName %in% c("unspInt", "behUnspInt")
    }
    )

    interactionEffects <-
      lapply(myeffects, function(x) {
             x[x$shortName %in% c("unspInt", "behUnspInt"), ]
    }
    )
    interactionEffectsl <-
      lapply(myeffects, function(x) {
             x$shortName %in% c("unspInt", "behUnspInt")
    }
    )
    ## store effects objects as we may need to recreate them
    f$interactionEffects <- interactionEffects
    f$basicEffects <- basicEffects
    f$interactionEffectsl <- interactionEffectsl
    f$basicEffectsl <- basicEffectsl
  # } else {
  #   myCompleteEffects <- ff$myCompleteEffects
  #   basicEffects <- ff$basicEffects
  #   interactionEffects <- ff$interactionEffects
  #   basicEffectsl <- ff$basicEffectsl
  #   interactionEffectsl <- ff$interactionEffectsl
  #   types <- ff$types
  # }
  ans <- .Call(C_effects, PACKAGE=pkgname, pData, basicEffects)
  pModel <- ans[[1]][[1]]
  for (i in seq(along=(ans[[2]]))) ## ans[[2]] is a list of lists of
    ## pointers to effects. Each list corresponds to one
    ## dependent variable
  {
    effectPtr <- ans[[2]][[i]]
    basicEffects[[i]]$effectPtr <- effectPtr
    interactionEffects[[i]]$effect1 <-
      basicEffects[[i]]$effectPtr[match(interactionEffects[[i]]$effect1,
                                        basicEffects[[i]]$effectNumber)]
        interactionEffects[[i]]$effect2 <-
          basicEffects[[i]]$effectPtr[match(interactionEffects[[i]]$effect2,
                                            basicEffects[[i]]$effectNumber)]
        interactionEffects[[i]]$effect3 <-
          basicEffects[[i]]$effectPtr[match(interactionEffects[[i]]$effect3,
                                            basicEffects[[i]]$effectNumber)]
  }
  ans <- .Call(C_interactionEffects, PACKAGE=pkgname,
               pModel, interactionEffects)
  ## copy these pointers to the interaction effects and then insert in
  ## effects object in the same rows for later use
  for (i in 1:length(ans[[1]])) ## ans is a list of lists of
    ## pointers to effects. Each list corresponds to one
    ## dependent variable
  {
    if (nrow(interactionEffects[[i]]) > 0) {
      effectPtr <- ans[[1]][[i]]
      interactionEffects[[i]]$effectPtr <- effectPtr
    }
    myCompleteEffects[[i]][basicEffectsl[[i]], ] <- basicEffects[[i]]
    myCompleteEffects[[i]][interactionEffectsl[[i]],] <-
      interactionEffects[[i]]
    ##myeffects[[i]] <- myeffects[[i]][order(myeffects[[i]]$effectNumber),]
  }
  ## remove the effects only created as underlying effects
  ## for interaction effects. first store the original for use next time
  myeffects <- lapply(myCompleteEffects, function(x) {
                      x[x$requested, ]
               }
  )
  if (!initC) {
    ans <- matrix(0, 1, length(myeffects$type != 'gmm'))
    ## create a grid of periods with group names in case want to
    ## parallelize using this or to access chains easily
    groupPeriods <- attr(f, "groupPeriods")
    z$callGrid <- cbind(rep(1:nGroup, groupPeriods - 1),
                        as.vector(unlist(sapply(groupPeriods - 1,
                                                function(x) 1:x))))
    }

  ##store address of model
  f$pModel <- pModel
  ans <- reg.finalizer(f$pModel, clearModel, onexit = FALSE)
  if (x$MaxDegree == 0 || is.null(x$MaxDegree)) {
    MAXDEGREE <-  NULL
  } else {
    MAXDEGREE <- x$MaxDegree
    storage.mode(MAXDEGREE) <- "integer"
  }
  if (x$UniversalOffset == 0 || is.null(x$UniversalOffset)) {
    UNIVERSALOFFSET <-  NULL
  } else {
    UNIVERSALOFFSET <- x$UniversalOffset
    storage.mode(UNIVERSALOFFSET) <- "double"
  }
  if (x$modelType == 0 || is.null(x$modelType)) {
    MODELTYPE <-  NULL
  } else {
    MODELTYPE <- x$modelType
    storage.mode(MODELTYPE) <- "integer"
  }
  if (x$behModelType == 0 || is.null(x$behModelType)) {
    BEHMODELTYPE <-  NULL
  } else {
    BEHMODELTYPE <- x$behModelType
    storage.mode(BEHMODELTYPE) <- "integer"
  }
  if (z$cconditional) {
    CONDVAR <- z$condname
    CONDTARGET <- attr(f, "change")
    ##   cat(CONDTARGET, "\n")
  } else {
    CONDVAR <- NULL
    CONDTARGET <- NULL
  }
  simpleRates <- TRUE
  if (any(!z$effects$basicRate & z$effects$type =="rate")) {
    simpleRates <- FALSE
  }
  z$simpleRates <- simpleRates
  ans <- .Call(C_setupModelOptions, PACKAGE=pkgname,
               pData, pModel, MAXDEGREE, UNIVERSALOFFSET, CONDVAR, CONDTARGET,
               profileData, z$parallelTesting, MODELTYPE, BEHMODELTYPE, z$simpleRates,
               x$normSetRates)
  if (x$maxlike) {
    # if (!initC) {
    #   ## set up chains and do initial steps
      types <- attr(f, "types")
      nbrNonMissNet <- attr(f, "numberNonMissingNetwork")
      nbrMissNet <- attr(f, "numberMissingNetwork")
      nbrNonMissBeh <- attr(f, "numberNonMissingBehavior")
      nbrMissBeh <- attr(f, "numberMissingBehavior")

      if (sum(nbrMissNet + nbrNonMissNet) > 0) {
        z$prmin <- nbrMissNet/ (nbrMissNet + nbrNonMissNet)
      } else {
        z$prmin <- rep(0, length(nbrMissNet))
      }
      if (sum(nbrMissBeh + nbrNonMissBeh) > 0) {
        z$prmib <-   nbrMissBeh/ (nbrMissBeh + nbrNonMissBeh)
      } else {
        z$prmib <- rep(0, length(nbrMissBeh))
      }
      ## localML
      if (is.null(x$localML)) {
        z$localML <- FALSE
      } else {
        z$localML <- x$localML
      }
      local <- ifelse(is.na(effects$local[effects$include]),
                      FALSE, effects$local[effects$include])
      if (z$localML & any(!local)) {
        stop("Non-local effect chosen.")
      }
      z$probs <- c(x$pridg, x$prcdg, x$prper, x$pripr, x$prdpr, x$prirms, x$prdrms)
      ans <- .Call(C_mlMakeChains, PACKAGE=pkgname, pData, pModel,
                   z$probs, z$prmin, z$prmib, x$minimumPermutationLength,
                   x$maximumPermutationLength, x$initialPermutationLength,
                   z$localML)
      f$minimalChain <- ans[[1]]
      f$chain <- ans[[2]]
    # } else ## set up the initial chains in the sub processes
    # {
    #   ans <- .Call("mlInitializeSubProcesses",
    #                PACKAGE=pkgname, pData, pModel,
    #                z$probs, z$prmin, z$prmib,
    #                x$minimumPermutationLength,
    #                x$maximumPermutationLength,
    #                x$initialPermutationLength, ff$chain,
    #                z$localML)
    #   f$chain <- ff$chain
    # }
  }
  f$simpleRates <- simpleRates
  f$myeffects <- myeffects
  f$myCompleteEffects <- myCompleteEffects
  # if (!initC) {
    # if (is.null(z$print) || z$print) DataReport(z, x, f)
    f$randomseed2 <- z$randomseed2
  # } else {
  #   f$randomseed2 <- ff$randomseed2
  # }
  f$observations <- attr(f, "observations") + 1
  f$depNames <- names(f[[1]]$depvars)
  f$groupNames <- names(f)[1:nGroup]
  f$nGroup <- nGroup
  f$types <- types
  if (!initC) {
    z$f <- f
    ##z <- initForAlgorithms(z)
    z$periodNos <- attr(data, "periodNos")
    z$f$myeffects <- NULL
    z$f$myCompleteEffects <- NULL
    if (!returnDeps) {
      z$f[1:nGroup] <- NULL
    }
  }
  if (initC || (z$int == 1 && z$int2 == 1 &&
                (is.null(z$nbrNodes) || z$nbrNodes == 1))) {
    f[1:nGroup] <- NULL
  }
  FRANstore(f) ## store f in FRANstore
  z$returnDeps <- returnDeps
  z$returnDepsStored <- returnDeps
  z$observations <- f$observations
  z$returnChains <- returnChains
  z$byGroup <- byGroup
  z$byWave <- byWave
  z$returnDataFrame <- returnDataFrame
  z$nDependentVariables <- length(z$f$depNames)
  # if (initC) {
  #   NULL
  # } else {
    z
  # }
}
