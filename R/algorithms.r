#******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snijders/siena/
# *
# * File: algorithms.r
# *
# * Description: Use to call simstats outside the Robbins-Monro framework.
# * Only works for one period at the moment, although the chains are returned
# * inside a group/period list structure for future expansion.
# *****************************************************************************/
#/**
##@algorithms algorithms use to call simstats outside R-M
algorithms <- function(data, effects, x, ...)
{
	## initialize
	z <- algorithmsInitialize(data, effects, x, ...)

	z$iter  <-  0

	finalLoop <- FALSE
	## library(lattice)

	z$thetasub <- z$thetasub + 1
	z$thetaHistory[z$thetasub, ] <- z$theta

	repeat ## done twice once for main body and once for final loop
	{
		repeat ## done numiter times plus final loop
		{
			########################################################
			## initialize for a new set of samples
			########################################################
			z$iter <- z$iter + 1
			##print(z$theta)

			if (z$useOptim)
			{
				if (!finalLoop)
				{
					z$nIter = z$optimSchedule[z$iter]
					## cat(z$iter, z$nIter, '\n')
				}
				z$importanceSamplingWeights <- rep(1, z$nIter)
			}
			cat(z$iter, z$nIter, '\n')
			#########################################################
			## get the next set of samples
			#########################################################
			if (finalLoop && z$variedThetas)
			{
				thetas <- data.matrix(z$thetaHistory)[1:z$thetasub, , drop=FALSE]
				if (z$thetasub > 20)
				{
					thetas <- thetas[-c(1:10), , drop=FALSE]
				}
				else
				{
					if (z$thetasub > 3)
					{
						thetas <- thetas[-c(1:2), , drop=FALSE]
					}
				}
				z <- getSamples(z, x, z$nIter, thetas=thetas)
			}
			else
			{
				z <- getSamples(z, x, z$nIter)
			}
			predictions <- colSums(colMeans(z$fra)) - z$targets

			iiter <- 0
			#######################################################
			## do update steps with these samples
			######################################################
			z$maxVar <- var(c(rep(0, z$nIter - 1), 1)) / 10
			repeat
			{
				iiter <- iiter + 1
				cat('iiter', iiter,'\n')
				if (z$useOptim)
				{
					## store old theta
					z$oldTheta <- z$theta

					z <- doOptimStep(z, x)
					if (z$thetasub < nrow(z$thetaHistory))
					{
						z$thetasub <- z$thetasub + 1
						z$thetaHistory[z$thetasub, ] <- z$theta
					}

					tmp <- getPredictions(z$theta, z, predictions=FALSE)
					z$importanceSamplingWeights <- tmp$ps
					varPs <-  tmp$varPs

					if (iiter > z$maxiiter || varPs > z$maxVar ||
						z$optimFn == 2)
					{
						break
					}
					if (all(abs(z$oldTheta - z$theta) < 0.001))
					{
						break
					}
				}
				else if (z$responseSurface)
				{
					z <- responseSurfaceChangeStep(z)
					print(z$theta)
					if (z$thetasub < nrow(z$thetaHistory))
					{
						z$thetasub <- z$thetasub + 1
						z$thetaHistory[z$thetasub, ] <- z$theta
					}
					break
				}
				else
				{
					z <- doAlgorithmChangeStep(z, x, predictions)
					print(z$theta)
					if (z$thetasub < nrow(z$thetaHistory))
					{
						z$thetasub <- z$thetasub + 1
						z$thetaHistory[z$thetasub, ] <- z$theta
					}
					cat(z$oldTheta - z$theta, '\n')

					if (finalLoop && all(abs(z$oldTheta - z$theta) < 0.01))
					{
						##browser()
						cat(z$oldTheta - z$theta, '\n')
						break
					}
					else
					{
						tmp <- getPredictions(z$theta, z)
						predictions <- tmp$predictions
						varPs <- tmp$varPs
						if (iiter > z$maxiiter || is.na(predictions))
						{
							break

						}
					}
				}
			}
			#######################################################
			## finished with these samples
			#######################################################
			if (z$thetasub > 1)
			{
				myformula <- as.formula(paste(z$formula, z$thetasub))
				print(xyplot(myformula, data=z$thetaHistory[1:z$thetasub, ],
						outer=TRUE,  scales="free",
						type="l"))
			}
			## clear the storage
			if (z$useHistory)
			{
				print(z$iter)
				z$storedLik0[[z$iter]] <- z$lik0
				if (z$iter > 3)
				{
					## do not remove an element but replace it by a null list
					z$storedLik0[z$iter - 3] <- list(NULL)
				}
				nbrStoredChains <- sum(sapply(z$storedLik0, length))
				print(sapply(z$storedLik0, length))
				keep <- nbrStoredChains + as.numeric(z$maxlike)
			}
			else
			{
				keep <- as.numeric(z$maxlike)
			}
			if (!z$useHistory || z$iter > 3)
			{
				if (z$nbrNodes > 1)
				{
					parSapply(z$cl, 1:nrow(z$callGrid), function(i, keep)
						{
							f <- FRANstore()
							.Call(C_clearStoredChains, PACKAGE=pkgname,
								f$pModel, keep, i)
						}, keep=keep)
				}
				else
				{
					sapply(1:nrow(z$callGrid), function(i)
						{
							f <- FRANstore()
							.Call(C_clearStoredChains, PACKAGE=pkgname,
								f$pModel, keep, i)
						})
				}
			}
			if (z$iter == z$numiter || finalLoop)
				## (iiter > z$maxiiter && varPs < z$maxVar / 2) || finalLoop )
			{
				##iter <- 0
				break
			}
		}  # end of inner repeat
		## only final loop to do now
		if (finalLoop)
		{
			break
		}
		z$diag <- FALSE
		z$gain <- z$gain / 2
		z$nIter <- z$finalIter
		z$numiter <- 1
		finalLoop <- TRUE
		if (z$useOptim)
		{
			z$maxiiter <- 0
		}
		else
		{
			z$maxiiter <- 20
		}
		z$maxit <- 100
		z$optimFn <- z$optimFinal
		z$useOptim <- z$useOptim && z$useOptimFinal
	}  ## end of outer repeat loop
	if (z$thetasub > 1)
	{
		myformula <- as.formula(paste(z$formula, z$thetasub))
		print(xyplot(myformula, data=z$thetaHistory[1:z$thetasub, ],
				outer=TRUE,  scales="free", type="l"))
	}

	z <- terminateFRAN(z, x)
	if (z$nbrNodes > 1)
	{
		stopCluster(z$cl)
	}
	if (z$useHistory)
	{
		z$storedLik0 <- NULL ## reduce size
	}
	z$zzz <- NULL
	z$FRAN <- NULL
	## if FRAN is there you will load RSiena on startup
	## if the return object is in your
	## workspace. Then it is difficult to recreate RSiena.
	z
}

##@doOptimStep algorithms Do one call to optim
doOptimStep <- function(z, x)
{
	theta <- z$theta
	theta[z$posj] <- log(theta[z$posj])
	if (z$optimFn == 1)
	{
		tmp <- optim(theta, optimFn, gr=derivFn, z=z,
			method='BFGS', control=list(maxit=z$maxit, trace=100))
	}
	else
	{
		tmp <- optim(theta, optimFn2, gr=derivFn2, z=z, method='BFGS',
			control=list(maxit=z$maxit, trace=100))
	}
	move <- sqrt(sum((theta - tmp$par)^2))
	maxmove <- z$gain
	if (move > maxmove)
	{
		move <- z$gain/move
	}
	else
	{
		move <- 1
	}
	z$par <- theta -  move * (theta - tmp$par)
	z$par[z$posj] <- exp(z$par[z$posj])
	z$theta <- z$par
	z
}

##@derivFn algorithms Gradient function for optim for sum(weighted(log ratios))
derivFn <- function(theta, z)
{
	theta[z$posj] <- exp(theta[z$posj])
	sc <- getLikelihoods(theta, z, getScores=TRUE)$sc
	sc <- rowMeans(sc)
	sc[z$posj] <- sc[z$posj] * theta[z$posj]
	if (z$useHistory && z$iter > 1) #was 5
	{
		start <- 1
		prevScores <-
			sapply(z$storedLik0, function(x)
				{
					if (is.null(x))
					{
						rep(0, length(z$theta))
					}
					else
					{
						## calculate where they were
						iterSequence <- start : (start + length(x) - 1)
						prevSc <- t(getLikelihoods(theta, z,
								getScores=TRUE,
								iterSequence= iterSequence)$sc)
						prevSc <- colMeans(prevSc)
						prevSc [z$posj] <- prevSc[z$posj] * theta[z$posj]
						start <<- start + length(x)
						prevSc
					}
				}
				)
		use <- max(1, z$iter - 3) : (z$iter - 1)
		use1 <- 1 : min(3, z$iter - 1)
		useWeights <- z$optimWeights[use1]
		sumWeights <- sum(useWeights)
		weights <- useWeights * z$optimWeight / sumWeights
		weight <- 1 - z$optimWeight
		weights <- rev(weights)
		sc <- weight * sc + rowSums(rep(weights, each=length(z$theta)) *
			prevScores[, use, drop=FALSE])
	}
	-sc
}

##@derivFn2 algorithms Gradient function for optim for sum(weighted(log ratios))
derivFn2 <- function(theta, z)
{
	opt <- optimFn2(theta,z)
	theta[z$posj] <- exp(theta[z$posj])
	tmp <- getLikelihoods(theta, z, getScores=TRUE)
	sc <- t(tmp$sc)
	sc[, z$posj] <- sc[, z$posj] * rep(theta[z$posj], each=nrow(sc))
	loglik <- tmp$lik
	ps <- exp(loglik - z$lik0)
	if (z$verbose)
	{
		stem(ps/sum(ps))
	}
	- 1/exp(opt) /z$nIter* colSums(sc*ps)
}
##@optimFn algorithms Function for optim for sum(weighted(log ratios))
optimFn <- function(theta, z)
{
	theta[z$posj] <- exp(theta[z$posj])
	loglik <- getLikelihoods(theta, z, getScores=FALSE)$lik
	loglik <- weighted.mean(loglik, z$importanceSamplingWeights)
	if (z$useHistory && z$iter > 1)
	{
		start <- 1
		prevLoglik <- sapply(z$storedLik0, function(x)
			{
				if (is.null(x))
				{
					0
				}
				else
				{
					## calculate where they were
					iterSequence <- 1 : ( start + length(x) - 1)
					tmp <- mean(getLikelihoods(theta, z,
							getScores=FALSE,
							iterSequence=
								iterSequence)$lik)
					start <<- start + length(x)
					tmp
				}
			}
			)
		use1 <- 1 : min(3, (z$iter - 1))
		use <- max(1, z$iter - 3) : (z$iter - 1)
		useWeights <- z$optimWeights[use1]
		sumWeights <- sum(useWeights)
		weights <- useWeights * z$optimWeight / sumWeights
		weight <- 1 - z$optimWeight
		weights <- rev(weights)
		#browser()
		loglik <- weight * loglik + sum(weights * prevLoglik[use])
	}
	- loglik
}
##@optimFn2 algorithms Function for optim for log(mean(likelihood ratios))
optimFn2 <- function(theta, z)
{
	theta[z$posj] <- exp(theta[z$posj])
	loglik <- getLikelihoods(theta, z, getScores=FALSE)$lik
	ps <- exp(loglik - z$lik0)
	-  log(mean(ps))
}

##@doGetProbabilitiesFromC algorithms Get likelihood for a stored chain
forwardGetProbabilitiesFromC <- function(index, z, getScores=FALSE)
{
	f <- FRANstore()
	anss <- apply(z$callGrid, 1, function(x)
		.Call(C_getChainProbabilities, PACKAGE = pkgname, f$pData,
			f$pModel, as.integer(x[1]), as.integer(x[2]),
			as.integer(index), f$myeffects, z$thetaMat[1,], getScores)
		)
	ans <- list()
	ans[[1]] <- sum(sapply(anss, "[[", 1))
	ans[[2]] <- sapply(anss, "[[", 2)
	ans[[3]] <- sapply(anss, "[[", 3)

	ans
}

##@getLikelihoods algorithms Get likelihoods from C for stored chainss
getLikelihoods <- function(theta, z, getScores=FALSE, iterSequence)
{
	z$thetaMat <- matrix(theta, nrow=1)
	if (z$maxlike || z$nbrNodes == 1)
	{
		if (missing(iterSequence))
		{
			iterSequence <- length(z$lik0) : 1
		}
		if (any(iterSequence < 0))
		{
			browser()
		}
		anss <-	lapply(iterSequence, function(i, z)
			getProbabilitiesFromC(z, i, getScores=getScores), z=z)
	}
	else
	{
		if (missing(iterSequence))
		{
			blocksize <- ceiling(z$nIter / z$nbrNodes)
			iterSequence <- blocksize : 1
			iterSequence <- rep(iterSequence, z$nbrNodes)[1:z$nIter]
		}
		anss <- parLapply(z$cl, iterSequence, forwardGetProbabilitiesFromC,
			z, getScores=getScores)
	}
	lik <- sapply(anss, "[[", 1)
	sc <- sapply(anss, "[[", 2)
	deriv <- sapply(anss, "[[", 3)
	list(lik=lik, sc=sc, deriv=deriv)

}

##@getPredictions algorithms Get predicted likelihood from stored simulations
getPredictions <- function(theta, z, predictions=TRUE)
{
	ps <- getLikelihoods(theta, z)$lik
	ps <- exp(ps - z$lik0)
	ps <- ps / sum(ps)
	if (z$verbose && any(!is.na(ps)))
	{
		stem(ps)
	}
	if (predictions)
	{
		if (is.na(var(ps)) || var(ps) > z$maxVar)
		{
			predictions <- NA
		}
		else
		{
			fras <- apply(z$fra, 3, function(x) colSums(x * ps))
			## dim fras is nwaves-1 by number of parameters (vector if nwaves=2)
			## so recreate it
			dim(fras) <- dim(z$fra)[2:3]
			predictions <- colSums(fras) - z$targets
		}
		list(predictions=predictions, varPs=var(ps), ps)
	}
	else
	{
		list(varPs=var(ps), ps=ps)
	}
}

##@algorithmsInitialize algorithms Initialize algorithm setup
algorithmsInitialize <-
	function(data, effects, x, scale=0.2, nIter=10,
		verbose=TRUE, numiter=10, maxiiter=10, useOptim=FALSE, diag=TRUE,
		finalIter=100, optimMaxit=1, optimWeight=0.7, useHistory=FALSE,
		optimSchedule=c(rep(c(20, 50, 100), c(10, 10, 10))),
		responseSurface=FALSE, variedThetas=FALSE,
		optimFinal=2, useOptimFinal=TRUE, nbrNodes=1,
		returnDataFrame=FALSE, clusterType=c("PSOCK", "FORK"))
	{
		##Report(openfiles=TRUE, type="n", silent=TRUE) #initialise with no file
		z  <-  NULL
		## get the function: simstats or maxlikec from RSiena or RSienaTest
		z$FRAN <- getFromNamespace(x$FRANname, pkgname)
		z$maxlike <- x$maxlike
		x$cconditional <-  FALSE
		z$gain <- scale
		z$diag <- diag
		z$nIter <- nIter
		z$finalIter <- finalIter
		z$verbose <- verbose
		z$useOptim <- useOptim
		z$useHistory <- useHistory
		z$maxiiter <- maxiiter
		z$numiter <-  numiter
		z$print <- FALSE
		z$responseSurface <- responseSurface
		z$maxit <- optimMaxit
		z$optimWeight <- optimWeight
		z$optimFinal <- optimFinal
		z$optimFn <- 1
		z$useOptimFinal <- useOptimFinal
		z$thetasub <- 0
		z$variedThetas <- variedThetas
		z$nbrNodes <- nbrNodes
		if (!is.null(x$randomSeed))
		{
			set.seed(x$randomSeed, kind="default")
		}
		else
		{
			if (exists(".Random.seed"))
			{
				rm(.Random.seed, pos=1)
				RNGkind(kind="default")
			}
		}


		z <- initializeFRAN(z, x, data, effects, prevAns=NULL, initC=FALSE,
			returnDeps=FALSE)

		if (z$nbrNodes > 1)
		{
			if (z$maxlike && z$observations < 3)
			{
				message("cannot use Multiple processors with ML with only 2 periods")
				z$nbrNodes <- 1
			}
			else
			{
				## require(parallel)
				clusterType <- match.arg(clusterType)
				if (clusterType == "PSOCK")
				{
					clusterString <- rep("localhost", nbrNodes)
					z$cl <- makeCluster(clusterString, type=clusterType,
						outfile = "cluster.out")
				}
				else
				{
					z$cl <- makeCluster(nbrNodes, type=clusterType,
						outfile="cluster.out")
				}
				clusterCall(z$cl, library, pkgname, character.only = TRUE)
				f <- FRANstore()
				clusterCall(z$cl, storeinFRANstore, f)
				clusterCall(z$cl, initializeFRAN, z, x, initC = TRUE,
					returnDeps=FALSE)
				clusterSetRNGStream(z$cl, iseed=as.integer(runif(1,
							max=.Machine$integer.max)))
				#	parLapply(z$cl, c('ff1','ff2'), sink)
			}
			if (z$maxlike)
			{
				z$int2 <- nbrNodes
			}
		}


		## more 'parameters' expected by simstats!
		z$Deriv <- TRUE
		if (z$nbrNodes == 1)
		{
			z$cl <- NULL
		}
		z$Phase <- 1

		z$addChainToStore <- TRUE
		if (useOptim)
		{
			z$thetaHistory <-
				matrix(NA, nrow=length(optimSchedule) * z$maxiiter + 100,
					ncol=length(z$theta))
		}
		else
		{
			z$thetaHistory <- matrix(NA, nrow=z$numiter*z$maxiiter + 100,
				ncol=length(z$theta))
		}
		colnames(z$thetaHistory) <- z$effects$shortName
		z$thetaHistory <- data.frame(z$thetaHistory)
		if (z$useOptim)
		{
			z$optimSchedule <- optimSchedule
			z$numiter <- length(z$optimSchedule)
			z$optimWeights <- z$optimWeight ^ (1:z$numiter)
			if (z$useHistory)
			{
				z$storedLik0 <- vector("list", z$numiter)
			}
		}
		z$formula <- paste(names(z$thetaHistory), collapse="+")
		z$formula <- paste(z$formula, "~ 1:")
		z
	}

##@doCreateChains algorithms Create storage for forward simulation chains
# doCreateChains <- function()
# {
# 	f <- FRANstore()
# 	ans <- .Call(C_createChainStorage, PACKAGE=pkgname,
# 		f$pData, f$pModel, f$simpleRates)
# 	f$pChain <- ans
# 	FRANstore(f)
# }

##@getSamples algorithms Get a set of samples from FRAN
getSamples <- function(z, x, nIter, thetas=NULL)
{
	## get a set of z$nIter samples
	sdf <- vector("list", nIter)
	if (is.null(thetas)) ## just use the theta on the object
	{
		doDeriv <- TRUE
		thetas <- matrix(z$theta, nrow=1)
	}
	else
	{
		doDeriv <- FALSE
	}
	if (z$nbrNodes > 1 && (!z$maxlike || nrow(thetas) > 1 || z$observations > 2))
	{
		zsmall <- getFromNamespace("makeZsmall", pkgname)(z)
		if (nrow(thetas) == 1) # straight repeats with same theta
		{
			if (z$maxlike)
			{
				## not a parLapply as we only parallelize by wave and
				## this is done in maxlikec automatically controlled by z$int2
				tmp <- lapply(1:nIter, function(i, z, FRAN)
					FRAN(z, NULL, returnLoglik=TRUE),
					z=zsmall, FRAN=z$FRAN)
				sdf <- lapply(tmp, function(x) x$dff)
				sc <- array(0, dim=c(nIter, z$observations - 1, length(z$theta)))
			}
			else
			{
				tmp <- parLapply(z$cl, 1:nIter, simstats0c, z=zsmall,
					returnLoglik=TRUE)
				sc <- t(sapply(tmp, function(x)x$sc))
				nobs <- dim(tmp[[1]]$fra)[1]
				dim(sc) <- c(nIter, nobs, z$pp)
			}
			fra <- t(sapply(tmp, function(x)x$fra))
			nobs <- dim(tmp[[1]]$fra)[1]
			dim(fra) <- c(nIter, nobs, z$pp)
			z$fra <- fra
			lik <- sapply(tmp, function(x)x$loglik)
			dim(lik) <- c(z$observations - 1, z$nIter)
			z$zzz <- lapply(tmp, function(x)x$chain)
			z$lik0 <- colSums(lik)
		}
		else ## multiple different thetas. Keep same thetas on one process
		{
			sc <- array(0, dim=c(nIter, z$observations - 1, length(z$theta)))
			niter <- ceiling(nIter / nrow(thetas))
			oldtheta <- z$theta
			##alternative: turn off z$int2 and use parRapply here.
			## need another function defined outside this one to reduce
			##traffic
			tmp <-
				apply(thetas, 1, function(th, z, n)
					{
						z$theta <- th
						tt <- lapply(1:n, function(i, z1)
							{
								maxlikec(z1, NULL, returnLoglik=TRUE)
							}, z1=z)
						tt
					}, z=zsmall, n=niter)

			## likelihood
			tmpa <- lapply(tmp, function(x) lapply(x, "[[", 12))
			tmpa <- do.call(c, tmpa)
			z$lik0 <- do.call(rbind, tmpa)
			z$lik0 <- rowSums(z$lik0)

			##fra
			fra <- lapply(tmp, function(i) lapply(i, function(x)x$fra))
			fra <- do.call(c, fra)
			fra <- array(unlist(fra),
				dim=c(z$observations - 1, z$pp, length(fra)))
			fra <- aperm(fra, c(3, 1, 2))
			##deriv
			sdf <- vector('list', niter)
			for (i in 1:nrow(thetas))
			{
				tmpa <- tmp[[i]]
				tmp1 <- lapply(tmpa, function(x)x$dff)

				sdf[((i-1) * niter + 1):(i * niter)] <- tmp1
			}

			## chains
			zzz <- lapply(tmp, function(y)
				lapply(y, function(x)x$chain))
			zzz <- do.call(c, zzz)
			z$theta <- oldtheta
		}
	}
	else ## just one process
	{
		zzz <- vector("list", nIter)
		fra <- array(0, dim=c(nIter, z$observations - 1, length(z$theta)))
		sc <- array(0, dim=c(nIter, z$observations - 1, length(z$theta)))
		z$startTheta <- z$theta
		thetasub <- 1
		lik0 <- matrix(NA, ncol=nIter, nrow=z$observations - 1)
		oldtheta <- z$theta
		for (i in (1:nIter))
		{
			## extract theta and set up index for next one
			z$theta <- thetas[thetasub, ]
			if (thetasub < nrow(thetas))
			{
				thetasub <- thetasub + 1
			}
			else
			{
				thetasub <- 1
			}
			## ##############
			## do iteration
			## ##############
			returnLoglik <- TRUE
			ans <- myfran(z, x, returnLoglik=returnLoglik)
			fra[i, , ] <- ans$fra ## statistics
			if (x$maxlike)
			{
				sdf[[i]] <- ans$dff
			}
			else
			{
				sc[i, , ] <- ans$sc
			}
			zzz[[i]] <- ans$chain
			##	browser()
			lik0[, i] <- ans$loglik
			##	lik0[, i] <- Reduce("+", Reduce("+", ans$loglik))
		}
		z$lik0 <- colSums(lik0)
		z$theta <- oldtheta
		z$zzz <- zzz
	}
	z$fra <- fra
	z$sc <- sc
	z$sd <- apply(apply(z$fra, c(1,3), sum), 2, sd)
	z$sdf <- sdf
	if (doDeriv)
	{
		if (!z$maxlike)
		{
			z$dfra <- getFromNamespace("derivativeFromScoresAndDeviations",
				pkgname)(z$sc, z$fra, , , , TRUE, )
		}
		else
		{
			z$dfra <- t(as.matrix(Reduce("+", z$sdf) / length(z$sdf)))
		}
		if (inherits(try(z$dinv <- solve(z$dfra), silent=TRUE), 'try-error'))
		{
			z$dinv <- NULL
		}
	}
	z
}

##@myfran models wrapper so profiler will include this separately
myfran <- function(z, x, ...)
{
	z$FRAN(z, x, ...)
}

##@doAlgorithmChangeStep algorithms Change step for algorithms
doAlgorithmChangeStep <- function(z, x, fra)
{
	z$oldTheta <- z$theta
	x$diag <- z$diag
	z <- doChangeStep(z, x, fra)
	z
}


##@responseSurfaceChangeStep algorithms Change step for linear surface
responseSurfaceChangeStep <- function(z, eps=0.1)
{
	npar <- length(z$theta)

	## make a grid
	mygrid <- matrix(1, ncol=npar, nrow=npar * 4)
	mygrid[cbind(seq(1, (4 * npar), 4), 1:npar)] <- 0.95
	mygrid[cbind(seq(2, (4 * npar), 4), 1:npar)] <- 1.05
	mygrid[cbind(seq(3, (4 * npar), 4), 1:npar)] <- 0.90
	mygrid[cbind(seq(4, (4 * npar), 4), 1:npar)] <- 1.10

	thetaForGrid <- z$theta

	thiseps <- ifelse(abs(z$theta) < 1e-10, eps, 0)
	myvals <- mygrid * rep(thetaForGrid + thiseps, each=nrow(mygrid))
	colnames(myvals) <- paste("theta", 1:length(z$theta), sep=".")
	## oldmaxVar <- z$maxVar
	z$maxVar <- 1
	predictedStats <- t(apply(myvals, 1, function(i, z)
			getPredictions(i, z)$predictions, z=z))
	colnames(predictedStats) <- paste("prediction", 1:length(z$theta), sep=".")
	mydf<- data.frame(myvals, pred=predictedStats, weights=rep(1, nrow(mygrid)))
	if (z$iter > 1)
	{
		mydf <- rbind(z$df2, mydf)
	}
	respCols <- (npar+1) : (2*npar)
	mylm <- lm(as.matrix(mydf[, respCols]) ~ . -weights ,
		data=mydf[, -respCols], weights=weights)
	coefs <- coef(mylm)[-1, ]
	newvals <- solve(t(coefs), - mylm$coef[1, ])
	move <- sqrt(sum((z$theta - newvals)^2))
	cat(move, diag(z$dfra), '\n')
	maxmove <- z$gain
	if (move > maxmove)
	{
		move <- z$gain/move
	}
	else
	{
		move <- 1
	}
	mydf$weights <- mydf$weights * z$optimWeight
	z$df2 <- mydf
	z$theta <- z$theta -  move * (z$theta - newvals)
	z
}

##@responseSurfaceQuadChangeStep algorithms Change step for quadratic surface
responseSurfaceQuadChangeStep <- function(z, eps=0.1)
{
	npar <- length(z$theta)

	## make a grid
	mygrid <- matrix(0, ncol=npar, nrow=npar * npar * 2)
	diag(mygrid[1:npar, ]) <- 1
	diag(mygrid[npar+1:npar, ]) <- -1
	mysubs <- matrix(0, nrow=(nrow(mygrid) - 2 * npar) * 2, ncol=2)
	mysubs[, 1] <- rep((2*npar + 1):nrow(mygrid), each=2)
	myrep <- rep(rep(1:(npar - 1), (npar - 1) : 1), 2)
	mysubs[, 2] <- t(cbind(myrep, myrep + rep(sequence((npar - 1):1), 4)))
	myrep <- c(rep(c(1.2, 0.8), each=nrow(mygrid)/ 2 - npar),
		rep(c(0.8, 1.2), (nrow(mygrid)/2 - npar)/2),
		rep(c(1.2, 0.8), (nrow(mygrid)/2 - npar)/2))
	myrep <- c(rep(c(1, -1), each=nrow(mygrid)/ 2 - npar),
		rep(c(-1, 1), (nrow(mygrid)/2 - npar)/2),
		rep(c(1, -1), (nrow(mygrid)/2 - npar)/2))
	mygrid[mysubs] <- myrep
	for (i in 1:npar)
	{
		ztheta <- z$theta
		mygrid[,i] <- mygrid[,i] * max(.1, (ztheta[i]+ eps)/10)
	}
	thetaForGrid <- z$theta
	#myvals <- mygrid * rep(thetaForGrid + eps, each=nrow(mygrid))
	myvals <-  mygrid + rep(thetaForGrid, each=nrow(mygrid))
	colnames(myvals) <- paste("theta", 1:length(z$theta), sep=".")
	##oldmaxVar <- z$maxVar
	z$maxVar <- 1
	predictedStats <- t(apply(myvals, 1, function(i, z)
			getPredictions(i, z)$predictions, z=z))
	names(predictedStats) <- paste("prediction", 1:length(z$theta), sep=".")
	ssdev <- apply(predictedStats, 1, function(x) sum(x^2))
	mydf<- data.frame(myvals, pred=ssdev/1000, weights=rep(1, nrow(mygrid)))
	if (z$iter > 1)
	{
		mydf <- rbind(z$df2, mydf)
	}
	myxx <- jitter(data.matrix(mydf[, 1:npar]))
	myxxx <- poly(myxx, degree=2, raw=TRUE)
	myxmat <- as.matrix(myxxx)
	mydf <- cbind(mydf, xx=myxmat)
	names(mydf)[(npar+3):ncol(mydf)] <-
		paste('x', 1:(ncol(mydf)-npar -2),sep='')
	mylm <- lm(pred ~ . -weights, data=mydf[,-(1:npar)], weights=mydf$weights)
	coefs <- coef(mylm)[-1]
	degrees <- attributes(myxxx)$degree
	## split the coefficients into A and b
	b <- coefs[degrees==1]
	AA <- coefs[degrees==2]
	A <- matrix(0, npar, npar)
	A[upper.tri(A, diag=TRUE)] <- AA
	A <- A + t(A)
	diag(A) <- diag(A)/2
	eigend <- eigen(A)
	if (any(eigend$values < 0))
	{
		mat1 <- eigend$vectors
		mat2 <- diag(eigend$values)
		mat2[mat2 < 0] <- 0
		aaa <- mat1 %*% mat2 %*% t(mat1)
		## require(MASS)
		newvals <- -ginv(aaa) %*% b/2
	}
	else
	{
		Ainv <- solve(A)
		newvals <- -Ainv %*% b/2
	}
	newvals <- newvals[, 1]
	move <- sqrt(sum((z$theta - newvals)^2))
	maxmove <- z$gain
	if (move > maxmove)
	{
		move <- z$gain/move
	}
	else
	{
		move <- 1
	}
	mydf$weights <- mydf$weights * z$optimWeight
	z$df2 <- mydf[, 1:(2*npar)]
	z$theta <- z$theta -  move * (z$theta - newvals)
	cat(move, newvals, z$theta,'\n')
	z
}

##@profileLikelihoods algorithm Calculate profile likelihood
profileLikelihoods <- function(resp, x, data, effects,
	i, j=NULL, gridl=c(0.8, 1.2), seqlen=5,
	maxit=2, method="BFGS", trace=0,
	nIter=100, ...)
{
	## initialize
	z <- algorithmsInitialize(data, effects, x, nIter=nIter,...)
	z$theta <- resp$theta
	z <- getSamples(z, x, nIter)
	z$importanceSamplingWeights <- rep(1, nIter)
	theta <- z$theta
	theta[z$posj] <- log(theta[z$posj])
	thetaFix <- theta
	fix <- rep(FALSE, length(theta))
	grid1 <- sort(seq(thetaFix[i]*gridl[1], thetaFix[i] * gridl[2], len=seqlen))
	fix[i] <- TRUE
	if (is.null(j)) # one dimensional density
	{
		zz <- sapply(grid1, function(x)
			{
				thetaFix[fix] <- x
				tmp <- optim(theta, profOptimFn, gr=profDerivFn,
					z=z, fix=fix, thetaFix=thetaFix,
					method=method,
					control=list(maxit=maxit, trace=trace))
				tmp$value
			}
			)
		if (z$posj[i])
		{
			grid1 <- exp(grid1)
		}
		plot(zz ~ grid1, type='b')
	}
	else # contour for 2
	{
		grid2 <- sort(seq(thetaFix[j]*gridl[1], thetaFix[j] * gridl[2],
				len=seqlen))
		fix[j] <- TRUE
		xy <- expand.grid(grid1, grid2)
		zz <- apply(xy, 1, function(x)
			{
				thetaFix[fix] <- x
				tmp <- optim(theta, profOptimFn, gr=profDerivFn,
					z=z, fix=fix, thetaFix=thetaFix,
					method=method,
					control=list(maxit=maxit, trace=trace))
				tmp$value
			}
			)
		zmat <- matrix(zz, nrow=length(grid1))
		if (z$posj[i])
		{
			grid1 <- exp(grid1)
		}
		if (z$posj[j])
		{
			grid2 <- exp(grid2)
		}
		contour(grid1, grid2, zmat)
		points(grid1, grid2)
	}
	z <- terminateFRAN(z, x)
	if (z$nbrNodes > 1)
	{
		stopCluster(z$cl)
	}
	z$zz <- zz
	invisible(z)
}

##@profOptimFn algorithms Function for optimizing profile likelihoods
profOptimFn <- function(theta, z, fix, thetaFix)
{
	theta[fix] <- thetaFix[fix]
	optimFn(theta, z)
}

##@profDerivFn algorithms Gradient function for optimizing profile likelihoods
profDerivFn <- function(theta, z, fix, thetaFix)
{
	theta[fix] <- thetaFix[fix]
	derivFn(theta, z)
}

##@clearStoredChains algorithms Clears storage used for chains in C.
#clearStoredChains <- function()
#{
#	f <- FRANstore()
#	.Call(C_clearStoredChains, PACKAGE=pkgname, f$pModel)
#}

##@doChangeStep algorithms change step for use in algorithms NB may be out of sync with phase 2
doChangeStep <- function(z, x, fra)
{
	## limit change.  Reporting is delayed to end of phase.
	if (is.null(x$maxmaxrat)){x$maxmaxrat <- 10.0}
	if (x$diag)## !maxlike at present
	{
		maxrat<- max(ifelse(z$sd, abs(fra)/ z$sd, 1.0))
		if (maxrat > x$maxmaxrat)
		{
			maxrat <- x$maxmaxrat / maxrat
			z$truncated[z$nit] <- TRUE
		}
		else
			maxrat <- 1.0
		fchange<- z$gain * fra * maxrat / diag(z$dfra)
	}
	else
	{
		fchange <- as.vector(z$gain * fra %*% z$dinv)
	}
	fchange[z$fixed] <- 0.0
	## check positivity restriction
	z$positivized[fchange > z$theta] <- z$positivized[fchange > z$theta] +1
	z$positivized[!z$posj] <- 0
	fchange <- ifelse(z$posj & (fchange > z$theta), z$theta * 0.5, fchange)

	z$theta <- z$theta - fchange
	if (!is.null(z$thav))
	{
		z$thav <- z$thav + z$theta
		z$thavn <- z$thavn + 1
	}
	z
}
