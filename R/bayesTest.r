##/*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: http://www.stats.ox.ac.uk/~snijders/siena
## *
## * File: bayesTest.r
## *
## * Description: This file contains the code to test results of sienaBayes,
## * and print results of the test.
## *
## ****************************************************************************/


##@simpleBayesTest Tests single parameters of sienaBayesFit objects
simpleBayesTest <- function(z, nfirst=z$nwarm+1, tested0=0,
	probs = c(0.025,0.975), ndigits=4){
	if (length(tested0) > 1)
	{
		stop('tested0 should be a number')
	}
	if (length(probs) != 2)
	{
		stop('probs should be a vector of length 2')
	}
	theEffects <- z$requestedEffects
	efNames <- format(paste(ifelse(theEffects$type == "creation",
				"creat", theEffects$type),
			theEffects$effectName)[!z$basicRate])
	if (all(theEffects$type[!z$basicRate] == 'eval'))
	{
		efNames <- format(theEffects$effectName[!z$basicRate])
	}
	else
	{
		efNames <- format(paste(ifelse(theEffects$type == "creation",
					"creat", theEffects$type),
				theEffects$effectName)[!z$basicRate])
	}
	credVal <- credValues(z, tested = tested0, theProbs = probs, nfirst=nfirst)
	mydf <- data.frame(matrix(NA, sum(!z$basicRate), 5))
	names(mydf) <- c(' ', 'varying', 'cred.from', '  cred.to', '  p  ')
	mydf[,1] <- efNames
	mydf[,2] <- ifelse(z$set2[!z$basicRate], '   -   ', '   +   ')
	mydf[z$fix[!z$basicRate],2] <- ' fixed '
	mydf[,3] <- format(round(credVal[!z$basicRate, 1], digits=ndigits))
	mydf[,4] <- format(round(credVal[!z$basicRate, 2], digits=ndigits))
	mydf[,5] <- format(round(credVal[!z$basicRate, 3], digits=4))
	mydf
}


##@multipleBayesTest Tests parameters of sienaBayesFit objects
multipleBayesTest <- function(z, testedPar, nfirst=z$nwarm+1, tested0=0, ndigits=4) {
	theEffects <- z$requestedEffects
	#	efNames <- format(paste(ifelse(theEffects$type == "creation",
	#									"creat", theEffects$type),
	#							theEffects$effectName)[!z$basicRate])
	if (all(theEffects$type[!z$basicRate] == 'eval'))
	{
		efNames <- format(theEffects$effectName[!z$basicRate])
	}
	else
	{
		efNames <- format(paste(ifelse(theEffects$type == "creation",
					"creat", theEffects$type),
				theEffects$effectName)[!z$basicRate])
	}
	# 7 lines borrowed from sienaBayes code:
	vec1 <- 4 - theEffects$randomEffects[theEffects$include]
	vec1[z$basicRate] <- 2
	vec1[z$ratePositions[[1]]] <- 1
	vec1[theEffects$fix[theEffects$include]] <- 5
	# now 1=rates group 1; 2=other rates; 3=randomly varying;
	# 4 = estimated non-varying (eta); 5 = non-estimated, fixed.
	# set1 = (1,2,3) set2 = (4)
	ntot <- sum(!is.na(z$ThinPosteriorMu[,1]))
	if (nfirst >= ntot-1)
	{
		stop('Warm sample too short')
	}
	nmax <- ntot-nfirst+1
	# Put mu and eta parameters in their proper places:
	p <- sum(vec1 %in% c(3,4,5))
	z$ThinObjective <- matrix(NA, nmax, p)
	vim <- vec1[vec1 %in% c(3,4,5)] == 3 # could be called z$muInObjective
	vie <- vec1[vec1 %in% c(3,4,5)] == 4 # could be called z$etaInObjective
	vif <- vec1[vec1 %in% c(3,4,5)] == 5 # could be called z$fixedInObjective
	z$ThinObjective[1:nmax, vim] <-
		z$ThinPosteriorMu[nfirst:ntot, z$objectiveInVarying, drop=FALSE]
	z$ThinObjective[1:nmax, vie] <- z$ThinPosteriorEta[nfirst:ntot,]
#	z$ThinObjective <- matrix(NA, nmax,
#		dim(z$ThinPosteriorEta)[2] +
#			dim(z$ThinPosteriorMu[,z$objectiveInVarying, drop=FALSE])[2])
#	vio <- vec1[vec1 %in% c(3,4)] == 3 # could be called z$varyingInObjective
#	z$ThinObjective[1:nmax, vio] <-
#		z$ThinPosteriorMu[nfirst:ntot, z$objectiveInVarying, drop=FALSE]
#	z$ThinObjective[1:nmax, !vio] <- z$ThinPosteriorEta[nfirst:ntot,]
#	p <- dim(z$ThinObjective)[2]
	if (inherits(testedPar, "matrix"))
	{
		A <- testedPar
		if (dim(A)[2] != p	)
		{
			stop(paste("If testedPar is a matrix, it should have",p,"columns."))
		}
		if (any(testedPar[,vif] != 0))
		{
			stop('Some fixed parameters were included in the tested set of effects.')
		}
		effeNames <- efNames
		testedNumber <- dim(A)[1]
		efNames <- rep(NA,testedNumber)
		for (i in (1:testedNumber))
		{
			efNames[i] <-
				paste(paste(A[i,],'*','parameter',1:p)[A[i,]!=0], collapse=" + ")
		}
		posteriorSample <- z$ThinObjective[,(vim|vie), drop=FALSE] %*% (t(A[,(vim|vie), drop=FALSE]))
		usedEffects <- apply(A,2,function(x){any(x!=0)})
		theUsedEffects <- paste(which(usedEffects),
			'. ',effeNames[usedEffects],sep='')
	}
	else
	{
		testedSet <- testedPar
		if ((min(testedSet) < 1) | (max(testedSet) > p))
		{
			stop(paste('testedPar should be a matrix or a set of numbers between 1 and', p))
		}
		testedNumber <- length(testedSet)
		efNames <- efNames[testedSet]
		if (any(vec1[testedSet] == 5))
		{
			stop('Some fixed parameters were included in the tested set of effects.')
		}
		posteriorSample <- z$ThinObjective[,testedSet, drop=FALSE]
		usedEffects <- (1:p) %in% testedPar
		theUsedEffects <- paste(which(usedEffects),'. ',efNames,sep='')
	}
	if (!(length(tested0) %in% c(1,testedNumber)))
	{
		stop(paste('tested0 must be one number or a vector of',
				testedNumber, 'numbers'))
	}
	if (length(tested0) == 1)
	{
		testedvec <- rep(tested0, testedNumber)
	}
	else
	{
		testedvec <- tested0
	}
	invCovParameters <- chol2inv(chol(cov(posteriorSample)))
	meanParameters <- colMeans(posteriorSample)
	standLength <- function(x){
		(x - meanParameters) %*% invCovParameters %*% (x - meanParameters)
	}
	posteriorStandLengths <- apply(posteriorSample, 1, standLength)
	boundary <- standLength(testedvec)
	postProb <- mean((posteriorStandLengths - rep(boundary, nmax)) > 0)
	result <- list(prob=postProb, chisquared=boundary,
		postDistances=posteriorStandLengths,
		nullValue=testedvec, effectNames=efNames, theUsedEffects=theUsedEffects,
		posteriorSample=posteriorSample)
	class(result) <- 'multipleBayesTest'
	result
}

print.multipleBayesTest <- function(x, descriptives=FALSE, ...){
	if (!inherits(x, 'multipleBayesTest'))
	{
		stop("not a legitimate multipleBayesTest object")
	}
	p <- length(x$effectNames)
	mydf <- as.data.frame(matrix(NA,p,4))
	if (p >= 2)
	{
		mydf[,1] <- paste('. ',x$effectNames,"= ")
	}
	else
	{
		mydf[,1] <- paste(x$effectNames,"= ")
	}
	mydf[,2] <- x$nullValue
	mydf[,3] <- round(colMeans(x$posteriorSample),4)
	mydf[,4] <- round(sqrt(diag(cov(x$posteriorSample))),4)
	names(mydf) <- c('', 'hyp', 'mean', 's.d')
	if (!(is.null(x$theUsedEffects)))
	{
		cat('Used effects:\n')
		for (i in 1:length(x$theUsedEffects)){cat(x$theUsedEffects[i],'\n')}
	}
	cat("\nTested hypothesis:\n")
	if (p >= 2) {print(mydf[,1:2])} else {cat(as.matrix(mydf[,1:2]),'\n')}
	if (descriptives)
	{
		cat('\nPosterior means and standard deviations:\n')
		if (p >= 2) {print(mydf[,c(1,3,4)])} else {cat(as.matrix(mydf[,c(1,3,4)]),'\n')}
		cat('\nPosterior covariance matrix:\n')
		print(round(cov(x$posteriorSample),6))
		cat('\n')
	}
	cat("\nposterior p-value:", round(x$prob, 4), "\n")
	cat("\ntest statistic:", round(x$chisquared, 2), "d.f. =",p,".\n\n")
	# Construct pattern of all sign combinations of p numbers
	s <- matrix(0,1,0)
	for (i in 1:p){
		pp <- 2^(i-1)
		s <- cbind(rbind(s,s), c(rep(1,pp),rep(-1,pp)))
	}
	propSignPattern <- function(si){
		mean(apply(x$posteriorSample, 1, function(x){identical(sign(x), si)}))
	}
	propSigns <- apply(s,1,propSignPattern)
	cat("Posterior proportion of sign patterns:\n")
	cat("", 1:p, "\n", sep="   ")
	pp <- dim(s)[1]
	# cbind transforms T to 1 and F to 0
	apply(cbind(s, propSigns), 1,
		function(y){cat(' ', ifelse((y[1:p] > 0.5), '>0 ', '<0 '),
			formatC(y[p+1], digits=3),'\n')})
	cat("\n")
	invisible(x)
}

stretch2 <- function(x){
	x[2] <- x[2] + 0.3*(x[2]-x[1])
	x
}


##@plot.multipleBayesTest methods
plot.multipleBayesTest <- function(x, xlim=NULL, ylim=NULL,	main=NULL, ...){
	if (!inherits(x, 'multipleBayesTest'))
	{
		stop("not a legitimate multipleBayesTest object")
	}
	# Makes a density plot of the posterior sample of distances.
	post <- x$postDistances
	par(oma=c(0,1,0,0), mar=c(5.1,5.1,4.1,2.1)) # to accommodate larger font for ylab
	# density by group
	if (is.null(xlim)){xlim <- stretch2(range(post))}
	xlim[1] <- min(0, xlim[1])
	xlim[2] <- max(xlim[2], x$chisquared + 0.1)
	# plot density so that truncation at 0 is clear;
	# since the distances are those from the posterior mean,
	# it is very unlikely that 0 is not in the support of the distribution.
	d1 <- density(c(post,-post))
	d1$x <- abs(d1$x)
	order.x <- order(d1$x)
	d1$x <- c(0,0,d1$x[order.x])
	d1$y <- c(0,d1$y[1],d1$y[order.x])
	if (is.null(ylim)){ylim <- stretch2(range(d1$y))}
	if (is.null(main)){main <- "posterior distances"}
	if (x$chisquared > xlim[2]){
		warning('note: observed chi-squared =',round(x$chisquared,1),
			'outside plot window.\n')}
	plot(d1, xlim=xlim, ylim=ylim, xlab='distance', ylab="density", col=4,
		cex=4, cex.lab=2, cex.main=2, main=main, lwd=2, ...)
	lines(c(x$chisquared, x$chisquared), c(0, max(d1$y)),lwd=2)
}

getNames <- function(x){
		# effect names without duplicated rate parameters
		# and with "unspecified interaction" replaced by
		# information about the effects in question
		b <- x$basicRate
		# True Parameters, i.e., all except rate parameters for groups 2 and up.
# earlier, can be deleted when checked:
#		tpar<- rep(NA,length(b))
#		for (i in (2:length(b))){tpar[i] <- !b[i]|(b[i]&!b[i-1])}
#		tpar[1] <- TRUE
		tpar <- !x$basicRate
		tpar[unlist(x$rateParameterPosition[[1]])] <- TRUE
		theNames <- x$requestedEffects$effectName[tpar]
		if (length(unlist(x$rateParameterPosition[[1]]) <= 1))
		{
		# Take away the ' (period 1)' in the first rate parameter
			theNames <- sub(' (period 1)','', theNames, fixed=TRUE)
		}
		theNames
}

##@extract.sienaBayes extracts samples from sienaBayesFit objects
extract.sienaBayes <- function(zlist, nfirst=zlist[[1]]$nwarm+1, extracted,
	sdLog=TRUE)
{
	if (!(is.list(zlist)))
	{
		stop('zlist must be a list of sienaBayesFit objects')
	}
	if (!all(sapply(zlist, function(z){inherits(z, "sienaBayesFit")})))
	{
		stop('all elements of zlist must be sienaBayesFit objects')
	}

	niter <- sapply(zlist,function(z){dim(z$ThinPosteriorMu)[1]})
	if (length(niter) > 1)
	{
		if (var(niter) > 1e-6)
		{
			stop("Length of these sienaBayes objects is not constant")
		}
	}
	nit <- niter[1] - nfirst + 1

	if (extracted %in% c("all", "rates")){
		# rate parameters
		EffName <- zlist[[1]]$effectName
		indices <- sort(which(!zlist[[1]]$generalParametersInGroup))
		nGroups <- zlist[[1]]$nGroup
		nind <- length(indices)
		res1 <- array(dim = c(nit, length(zlist), nGroups*length(indices)))
		fName1 <- rep('',dim(res1)[3])
		for (h in 1:length(zlist)){for (i in 1:nGroups){
			mini <- (i-1)*nind + 1
			maxi <- i*nind
			res1[,h,mini:maxi] <-
				zlist[[h]]$ThinParameters[(niter[h]-nit+1):niter[h],
			i,indices,drop=FALSE]
			fName1[mini:maxi] <- paste(EffName[indices],i)
		}}
	}

	if (extracted %in% c("all", "objective", "varying")){
		# mu
		varyings <- zlist[[1]]$varyingInEstimated
		nind <- sum(varyings)
		res2 <- array(dim = c(nit, length(zlist), (2*nind)))
		if (nind <= 0)
		{
			cat("Note: no varying parameters.\n")
			fName2 <- NULL
		}
		else
		{
			EffName <- getNames(zlist[[1]])[varyings]
			for (h in 1:length(zlist)){
				res2[,h,1:nind] <-
					zlist[[h]]$ThinPosteriorMu[(niter[h]-nit+1):niter[h], , drop=FALSE]
				postsig <- zlist[[h]]$ThinPosteriorSigma
				postsigg <- matrix(NA,dim(postsig)[1],nind)
				if (sdLog){
					for (k in 1:nind){postsigg[,k] <- 0.5*log(postsig[,k,k])}
				}
				else {
					for (k in 1:nind){postsigg[,k] <- sqrt(postsig[,k,k])}
				}
				res2[,h,(nind+1):(2*nind)] <-
					postsigg[(niter[h]-nit+1):niter[h],,drop=FALSE]
			}
			fName2 <- rep('',2*nind)
			fName2[1:nind] <- paste(EffName,'mu')
			if (sdLog){
				fName2[(nind+1):(2*nind)] <- paste(EffName,'log(s.d.)')
			}
			else {
				fName2[(nind+1):(2*nind)] <- paste(EffName,'s.d.')
			}
		}
	}

	if (extracted %in% c("all", "objective", "non-varying")){
		# eta
		nind <- sum(!zlist[[1]]$varyingParametersInGroup)
		res3 <- array(dim = c(nit, length(zlist), nind))
		if (nind <= 0)
		{
			cat("Note: no non-varying parameters.\n")
			fName3 <- NULL
		}
		else
		{
			EffName <- getNames(zlist[[1]])[!zlist[[1]]$varyingParametersInGroup]
			for (h in 1:length(zlist)){
				res3[,h,1:nind] <-
					zlist[[h]]$ThinPosteriorEta[(niter[h]-nit+1):niter[h], ,drop=FALSE]
			}
			fName3 <- EffName
		}
	}

	res <- NULL
	switch(extracted,
		'all' = {
			res <- array(c(res1,res2,res3),
				dim=c(dim(res1)[1], dim(res1)[2],
					dim(res1)[3]+dim(res2)[3]+dim(res3)[3]))
			fName <- c(fName1, fName2, fName3)
		},
		'rates' = {
			res <- res1
			fName <- fName1
		},
		'varying' = {
			res <- res2
			fName <- fName2
		},
		'non-varying' = {
			res <- res3
			fName <- fName3
		},
		'objective' = {
			res <- array(c(res2,res3),
				dim=c(dim(res2)[1], dim(res2)[2],
					dim(res2)[3]+dim(res3)[3]))
			fName <- c(fName2, fName3)
		},
		warning('wrong parameter given: extracted =',extracted,'\n')
	)
	if (!is.null(res)){
		dimnames(res) <- list(1:dim(res)[1], paste('chain', 1:length(zlist)),
			fName)
		if (any(is.na(res))){
			warning('warning: the extracted array contains NA elements\n')
		}
	}
	res
}


##@extract.posteriorMeans extracts posterior means from sienaBayesFit object
extract.posteriorMeans <- function(z, nfirst=z$nwarm+1, pmonly=1,
								excludeRates=FALSE, verbose=TRUE){
# produces a matrix with the groups in the rows
# and all effects in the columns, with for each effect
# first the posterior mean ("p.m.") and then the posterior standard deviation ("psd.")
	if (!inherits(z, "sienaBayesFit"))
	{
		stop('z must be a sienaBayesFit object')
	}
	ntot <- max(which(!is.na(z$ThinPosteriorMu[,1])))
	nit <- ntot - nfirst + 1
	if (excludeRates)
	{
		nind <- sum(z$objectiveInVarying)
	}
	else
	{
		nind <- sum(z$varyingParametersInGroup)
	}
	res <- matrix(NA, z$nGroup, 2*nind)
	if (nind <= 0)
	{
		cat("Note: no varying parameters. Null matrix is produced.\n")
	}
	else
	{
		EffName <- getNames(z)[z$varyingParametersInGroup]
		if (excludeRates)
		{
			EffName <- EffName[z$objectiveInVarying]
		}
		if (verbose)
		{
			cat(z$nGroup, ' groups\n')
		}
		for (h in 1:z$nGroup){
			if (verbose)
			{
				cat('.')
				flush.console()
				if ((h %% 50) == 0)
				{
				 cat(h, '\n')
				}
			}
			df <- sienaFitThetaTable(z, fromBayes=TRUE, tstat=FALSE,
									groupOnly=h, nfirst=nfirst)$mydf
			if (excludeRates)
			{
				seth <- which(z$varyingObjectiveParameters)
			}
			else
			{
				seth <- sort(union(z$ratePositions[[h]],
							which(z$varyingObjectiveParameters)))
			}
			posttheta <- df[seth,"value"]
			postsd <- df[seth,"se"]
			res[h,1:nind] <- posttheta
			res[h,(nind+1):(2*nind)] <- postsd
		}
		fName <- rep('',2*nind)
		fName[1:nind]  <- paste('p.m.',EffName)
		fName[nind + (1:nind)]  <- paste('psd.',EffName)
		dimnames(res) <- list(1:dim(res)[1], fName)
	}
	if (pmonly == 1)
	{
		res <- res[,1:nind]
	}
	else if (pmonly >= 2)
	{
		res <- res[,nind + (1:nind)]
	}
	if (verbose)
	{
		cat('*\n')
	}
	res
}


##@plotPostMeansMDS MDS plot of posterior means for sienaBayesFit object
plotPostMeansMDS <- function(x, pmonly=1, excludeRates=TRUE, nfirst=NULL, ...){
# This function makes an MDS plot of the posterior means in z;
# for the method: see MASS (book) p. 308.
# if pmonly=0 posterior means and standard deviations,
# if pmonly=1 only the posterior means,
# if pmonly=2 only the posterior standard deviations.
	if (!inherits(x, "sienaBayesFit"))
	{
		stop('x must be a sienaBayesFit object')
	}
	objectName <- deparse(substitute(x))
	if (is.null(nfirst))
	{
		nfirst <- x$nwarm+1
	}
#	requireNamespace(MASS)
	is.even <- function(k){k %% 2 == 0}
	is.odd <- function(k){k %% 2 != 0}
	message('extracting posterior means ...')
	pm <- extract.posteriorMeans(x, nfirst=nfirst, pmonly=pmonly,
									excludeRates=excludeRates)
	if (pmonly <= 0)
	{
		vars <- (1:dim(pm)[2])
	}
	else if (pmonly == 1)
	{
		vars <- is.odd(1:dim(pm)[2])
	}
	else
	{
		vars <- is.even(1:dim(pm)[2])
	}
	message('calculating MDS solution ...')
	corpm <- cor(t(pm[,vars]))
	mds <- isoMDS(1-corpm)
	plot(mds$points, type='n', main=paste('MDS',objectName),
											xlab='', ylab='', ...)
	text(mds$points, labels=as.character(1:dim(mds$points)[1]), ...)
	invisible(list(corpm=corpm, points=mds$points))
}