#' This file is the entry point for the "unit" tests. Currently there are only
#' before/after or persistence tests comparing values to previously stored
#' values. To record the current state run:
#'   R -e "require(RSienaTest); record_values <- T; RSienaTest:::run_tests(dir='.')"
#' to test against a saved state run:
#'   R -e "require(RSienaTest); RSienaTest:::run_tests(dir='.')"

#' Runs matching tests comparing the results to previously stores values
#' @param funcRegexp a regular expression for matching test function names.
#' @param dir directory name containing the test data (defaults to the
#' unitTests directory inside the package)
#' @return `T` iff all passed
run_tests <- function(funcRegexp=NA, dir=NA) {
	# change working directory (TODO: are all library paths accessible on all systems?)
	if (is.na(dir)) dir <- file.path(find.package('RSienaTest'), 'unitTests')
	setwd(dir)
	# select test by regular expression
	if (is.na(funcRegexp)) funcRegexp='^test.+'
	else funcRegexp=paste('^test.*', funcRegexp, '.*', sep='')
	# test suite
	suite <- RUnit::defineTestSuite('all',
		dirs=file.path(find.package(pkgname), 'unitTests'),
		testFileRegexp='.+test\\.R',
		testFuncRegexp=funcRegexp)
	res <- RUnit::runTestSuite(suite, verbose=T)$all$sourceFileResults
	# summary
	passed <- unlist(as.vector(sapply(names(res),
				function(file) {
					cat('\n', file, '\n', sep='')
					as.vector(sapply(names(res[[file]]),
							function(method) {
								test <- res[[file]][[method]]
								if (test$kind == 'success') {
									cat(' - ', method, '... passed\n', sep='')
									T
								} else {
									cat(' - ', method, '... ', test$kind, '\n', sep='')
									cat(test$msg, '\n', sep='')
									F
								}
							}))
				})))
	cat('passed: ', length(which(passed==T)),
		', failed: ', length(which(passed==F)), '\n', sep='')
	invisible(all(passed))
}

