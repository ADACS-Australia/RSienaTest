#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snijders/siena
# *
# * File: RSienaAutoHelpTable.r
# *
# * Description: This module contains the code for the automatic extraction of
# * documentation from the .Rd files
# *****************************************************************************/
##@RSienaAutoHelpTable Documentation Create latex source for a function table
RSienaAutoHelpTable <- function(RSienaDir="../RSiena/man")
{
	## Documentation of the Rd processing functions is in parseRd.pdf
	## a CRAN document.

	## All functions to be included need to be here in order, with
	## stage of processing relevant, and T or F in the arguments column
	## to indicate whether the arguments are to be included in the final
	## column of the table - TTT means all three arguments...
	RSienaRequiredFunctions <- data.frame(functionName=c(
			"sienaNodeSet",
			"sienaDependent",
			"coCovar",
			"varCovar",
			"coDyadCovar",
			"varDyadCovar",
			"sienaCompositionChange",
			"sienaDataCreate",
			"sienaDataCreateFromSession",
			"sienaGroupCreate",
			"effectsDocumentation",
			"getEffects",
			"print.sienaEffects",
			"includeEffects",
			"includeInteraction",
			"setEffect",
			"print01Report",
			"`",
			"sienaAlgorithmCreate",
			"siena07",
			"sienaFit",
			"sienaTimeTest"),
		stage=c("1", "3--5", rep("3", 10),
			rep("4", 6), "5", "5", "5", "6", "6"),
		arguments=c("", "F", "TTT",
			rep("F", 20)), stringsAsFactors=FALSE)

	numberPages <- nrow(RSienaRequiredFunctions)
	usage <- vector("list", numberPages)
	desc <- vector("list", numberPages)
	exam <- vector("list", numberPages)
	arguments <- vector("list", numberPages)

	## library(tools)

	for (i in 1:nrow(RSienaRequiredFunctions))
	{
		## find the path to the help file
		infile <- file.path(RSienaDir, paste(RSienaRequiredFunctions[i, 1],
				".Rd", sep=""))
		## parse it
		Rd <- tools::parse_Rd(infile)
		## extract the tags from the Rd object
		# used to be:		tags <- tools:::RdTags(Rd)
		# replaced by next two lines (to avoid :::)
		tags <- sapply(Rd, attr, "Rd_tag")
		if (!length(tags)) {tags <- character()}
		usage[[i]] <- Rd[[which(tags=="\\usage")]]
		desc[[i]] <- Rd[[which(tags=="\\description")]]
		exam[[i]] <- Rd[[which(tags=="\\examples")]]
		if (length(which(tags=="\\arguments")))
		{
			arguments[[i]] <- Rd[[which(tags=="\\arguments")]]
		}
	}

	## do some reformatting: use debug() to see what happens
	usage1 <- sapply(usage, function(x)paste(unlist(x), collapse=""))
	desc1 <- sapply(desc, function(x)paste(unlist(x), collapse=""))
	exam1 <- sapply(exam, function(x)paste(unlist(x), collapse=""))
	args1 <- lapply(arguments, function(x)
		{
			if (!is.null(x))
			{
				# used to be:                       tmp <- tools:::RdTags(x)
				# replaced by next two lines (to avoid :::)
				tmp <- sapply(x, attr, "Rd_tag")
				if (!length(tmp)) {tmp <- character()}
				tmp2 <- x[which(tmp=="\\item")]
				lapply(tmp2, function(x)list(x[[1]],
						gsub("\\n", "",
							x[[2]])))
			}
		})

	## optional reformatting: some examples which could be used

	##usage2 <- sub("getDocumentation=FALSE", "", usage1)
	##usage2 <- sub(", )", ")", usage2)
	##usage2 <- gsub("\\n", " ", usage2)
	##desc2 <- gsub("\\n", " ", desc1)
	##exam2 <- gsub("\\n", " ", exam1)
	##usage2 <- gsub("#", " ", usage2)
	##exam2 <- gsub("#", "\\#", exam2, fixed=TRUE)
	##usage2 <- gsub("$", "\\$", usage2, fixed=TRUE)

	##exam2 <- gsub("<", "\\textless", exam1)

	##usage2 <- sub("^ *", "", usage2)
	##exam2 <- sub("^ *", "", exam2)

	##@sortout Internal RSienaAutoHelpTable
	##more formatting
	sortout <- function(x, code=TRUE)
	{
		## remove trailing \n
		x <- gsub("\n$", "", x)
		## remove starting \n
		x <- gsub("^\n", "", x)
		## retain the others if after a )
		x <- gsub(")\n", ")\\\\ ", x)
		## remove the others
		x <- gsub("\n", " ", x)
		## escape hash
		x <- gsub("#", "\\#", x, fixed=TRUE)
		## escape dollar
		x <- gsub("$", "\\$", x, fixed=TRUE)
		##remove spaces
		x <- gsub("^ *", "", x)
		##remove spaces
		x <- gsub(" +", " ", x)
		if (code)
		{
			x <- paste("\\texttt{", x, "}", sep="")
		}
		x
	}
	## create the table: check doc.tex for the output
	mytab <-
		sapply(1:numberPages, function(i, rs, usage, exam, desc, arg)
			{
				x <- rs[i, ]
				includeArguments <- strsplit(x$arguments, "")[[1]]
				includeArguments <- as.logical(includeArguments)
				arg <- args1[[i]]
				arg <- arg[includeArguments]
				## sort out usage
				usage2 <- sortout(usage1[i])
				usage2 <- sub("getDocumentation=FALSE", "", usage2)
				if (nchar(x$functionName) > 20)
				{
					#browser()
					usage2 <- paste("\\\\ &&", usage2)
				}
				exam2 <- sortout(exam1[[i]])
				desc2 <- sortout(desc1[[i]], code=FALSE)
				arg <- sapply(arg, function(x)
					{
						paste('\\texttt{"', unlist(x[[1]]), '"}=',
							paste(unlist(x[[2]]), sep=" ", collapse=""),
							sep="", collapse="")
					})
				tmp <- paste(x$stage, "&", x$functionName, "&",
					usage2, "&",
					exam2, "&", desc2, paste(arg, sep=" ",
						collapse=""),
					"\\\\ \\hline",
					sep="", collapse="")
				print(tmp)
				tmp
			}, rs=RSienaRequiredFunctions, usage=usage, desc=desc,
			exam=exam, arg=args1)
	write.table(mytab, "docn.tex", row.names=FALSE, col.names=FALSE,
		quote=FALSE)
	## now run LaTeX on auto.tex in the doc directory.
}
##@RSienaSlowTest Tests Combine all examples on help pages and run them
RSienaSlowTest <- function(RSienaDir="../RSiena/", lib="library(RSiena)")
{
	helpdir <- dir(paste(RSienaDir, "/man", sep=""), pattern="Rd$",
		full.names=TRUE)
	tmp <- lapply(helpdir, function(x)
		{
			if (!grepl("maxlikefn", x)) ## obsolete, not-working function
			{
				tf <- tempfile("RSiena")
				tools::Rd2ex(x, tf)
				if (file.exists(tf))
				{
					on.exit(unlink(tf))
					tmpa <- readLines(tf)
					tmpb <- sub("##D", "", tmpa)
					tmpb <- ifelse (grepl("fix(", tmpb, fixed=TRUE),
						paste("#", tmpb), tmpb)
					if (basename(x) == "sienaDataCreateFromSession.Rd")
					{
						chnewdir <- paste("setwd(\"",
							RSienaDir,
							"/inst/examples\")", sep="")
						tmpb <- c("curdir <- getwd()",
							chnewdir, tmpb, "setwd(curdir)")
					}
					unlink(tf)
					c(paste("#", x), tmpb)
				}
			}
		}
		)
	tmp1 <- do.call(c, tmp)
	tmp1 <- c(lib, tmp1)
	tf <- tempfile("RSiena")

	write(tmp1, tf)
	source(tf, echo=TRUE)
	unlink(tf)
}
