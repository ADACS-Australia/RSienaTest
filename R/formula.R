#' Environment with effect definitions.
#'
#' Each effect is a function modifying the siena effects object.
#'
#' Effect aguments:
#' - eff: siena effects object
#' - var: dependent variable name
#'
#' Effect return value: The modified effects object.
local({
	# TODO: lazy loading, allEffects is not accessible right now
	# for (e in allEffects$shortName[allEffects$effectGroup %in% c('symmetricObjective')]) {
	# define a bunch of effects
	for (e in c('density', 'transTrip')) {
		assign(e, function(eff, var) {
			setEffect(eff, e, name=var, character=T)
		})
	}
	avAlt <- function(eff, var, ...) {
		arg <- as.list(sys.call())[-1]  # arguments as called
		unnamed <- arg[which(names(arg) == '')]  # unnamed arguments
		net <- as.character(unnamed[[1]]) # original object name
		setEffect(eff, avAlt, name=var, interaction1=net)
	}
}, envir=effectsenv <- new.env())

#' Constructs a siena model specification.
#'
#' @param ... list of formulas
#' @return object of class 'siena_formula'
siena_formula <- function(...) {
	# assert formula
	if (!all(sapply(m, has_lhs.formula))) stop('formula without response')
	# take response variable as name
	names(m) <- as.character(sapply(m, function(f) f[[2]]))
	structure(m, class='siena_formula')
}

print.siena_formula <- function(sformula) {
	for (f in sformula) {
		s <- as.character(f)
		cat(s[[2]], '~', s[-c(1, 2)], '\n')
	}
	invisible(f)
}

is.siena_formula <- function(f) {
	inherits(f, 'siena_formula')
}

#' Concatenates formulas for the same response variable.
#'
#' @param lhs object of class 'siena_formula'
#' @param rhs object of class 'siena_formula'
#' @return object of class 'siena_formula'
`+.siena_spec` <- function(lhs, rhs) {
	vars <- unique(c(names(lhs), names(rhs)))
	print(vars)
	vars <- sapply(vars, function(n) {
		print(lhs[[n]])
		print(rhs[[n]])
		if (is.null(lhs[[n]])) {
			return(rhs[[n]])
		} else if (is.null(rhs[[n]])) {
			return(lhs[[n]])
		} else {
			return(c.formula(lhs[[n]] + rhs[[n]]))
		}
	})
	do.call(siena_formula, vars)
}

# a <- siena_formula(x ~ density)
# b <- siena_formula(z ~ avAlt(x))

has_lhs.formula <- function(f) length(f) == 3

rhs.formula <- function(f) rhs <- f[[length(f)]]

lhs.formula <- function(f) if (has_lhs.formula(f)) return(f[[2]])

c.formula <- function(l, r) {
	if (lhs.formula(l) != lhs.formula(r)) stop('responses do not match')
	l[[length(l)]] <- call('+', rhs.formula(l), rhs.formula(r))
	l
}

#' Constructs a siena data object from a model specification.
#'
#' NOTE: this is seems all a bit hacky right now.  which argument is a data
#' object.  should be clearer which environment is used.  lots of error
#' checking to be added
#'
#' @param sformula a model specification
#' @return result of sienaDataCreate
extract_data <- function(sformula, allowOnly=F) {
	# assert class
	if (!is.siena_formula(sformula)) stop('not a model formula')
	# collect data objects
	dat <- list()
	for (f in sformula) {
		# add response variable
		response <- as.character(f[[2]])
		dat[[response]] <- get(response, envir=attr(f, '.Environment'))
		# other variables can appear as argument to calls
		vars <- as.list(attr(f, 'variables'))[-1]
		for (cal in vars) {
			if (is.call(cal)) {
				for (i in 2:length(cal)) {
					arg <- cal[[i]]
					dat[[as.character(arg)]] <- eval(arg)
				}
			}
		}
	}
	# convert common types to siena types
	for (n in names(dat)) {
		x <- dat[[n]]
		if (class(x) == 'array') {
			if (length(dim(x)) == 3) {
				dat[[n]] <- sienaDependent(x, allowOnly=allowOnly)
			}
		} else if (class(x) == 'matrix') {
			dat[[n]] <- sienaDependent(x, type='behavior', allowOnly=allowOnly)
		}
	}
	do.call(sienaDataCreate, dat)
}

#' Runs the effect definitions.
#'
#' NOTE: think about the arguments of effects, named/unnamed, quoted/raw
#'
#' @param sformula a model specification
#' @param eff siena effects object
#' @return siena effects object
run_effects <- function(sformula, eff) {
	if (!is.siena_formula(sformula)) stop('not a model formula')
	#
	for (f in sformula) {
		factors <- attr(f, 'factors')
		vars <- attr(f, 'variables')
		response <- as.character(f[[2]])
		for (col in 1:ncol(factors)) {
			# index of variables in this factor
			factor_col <- 1 + which(factors[,col] == 1)
			if (length(factor_col) == 1) {
				active <- vars[[factor_col]]
				# convert shorthand symbols to call form
				if (is.symbol(active)) active <- as.call(list(active))
				# add mandatory dependent variable and effects object
				active$var <- response
				active$eff <- eff
				# evaluate
				eff <- eval(active, envir=effectsenv)
			} else {
				stop('interaction not yet implemented')
			}
		}
	}
	eff
}

#' Convenience wrapper to siena07.
#'
#' @param sformula a model specification as produced by `siena_formula`
#' @param control algorithm specification as produced by `model.create`
#' @return result of the siena07 call
siena17 <- function(sformula, control, dat=NULL) {
	if (!is.siena_formula(sformula)) stop('not a model formula')

	# auto detect data from the formulas
	if (is.null(dat)) {
		dat <- extract_data(sformula)
	}

	eff <- getEffects(dat)
	eff <- run_effects(sformula, eff)
	siena07(control, data=dat, effects=eff, batch=TRUE)
}
