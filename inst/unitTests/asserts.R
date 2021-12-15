# This file contains test functions for comparing simulation outcomes commonly
# use by tests.
if (!exists('record_values') || !record_values) {
  record_values <- F
  cat('test file sourced in check mode!\n')
} else {
  cat('test file sourced in recording mode!\n')
}
cat('test directory:', getwd(), '\n')

#' @param name a name identifying the test
#' @return a valid file name with .Rd ending
file_name <- function(name) paste(gsub('[ .]', '_', name), '.Rd', sep='')

skip_recording <- function(name) {
  if (record_values && file.exists(file_name(name))) {
    cat('skip recording\n')
    T
  } else {
    F
  }
}

#' checks targets and parameter values of an estimated model
#' @param name a unique name for the value. most convenient if only one value
#' per function, then use  `match.call()[[1]]`.
#' @param value an estimated siena07 model
check_model_persistence <- function(name, value) {
  # save/load actual/expected values
  if (record_values) {
    save(value, file=file_name(name), ascii=T)
    return() # return is not a keyword but a function!!!
  }
  actual <- value
  load(file=file_name(name))
  expected <- value
  # first line break (runit does not append it after the '...')
  cat('\n')
  RUnit::checkEquals(expected$n1, actual$n1)
  # check model size (parameters and statistics (the latter for the gmm))
  RUnit::checkEquals(length(expected$theta), length(actual$theta))
  RUnit::checkEquals(length(expected$targets), length(actual$targets))
  if (is.null(expected$sienafit)) {
    RUnit::checkEquals(length(expected$effects), length(actual$effects))
    # permutation of actual effects to expected effects
    atoe <- sapply(expected$effects$shortName,
                   function(n) which(n == actual$effects$shortName)[[1]])
    # TODO: duplicate effect names (rates, same effects on different vars)
    # print(atoe)
    print(data.frame(expected.name=expected$effects$shortName,
                     # actual.name=actual$effects$shortName[atoe],
                     expected.targets=expected$targets,
                     diff.targets=abs(expected$targets - actual$targets[atoe]),
                     expected.theta=expected$theta,
                     diff.theta=abs(expected$theta - actual$theta[atoe])))
    RUnit::checkEquals(expected$effects$shortName, actual$effects$shortName[atoe])
    RUnit::checkEquals(expected$targets, actual$targets[atoe], tolerance=1e-12)
    RUnit::checkEquals(expected$theta, actual$theta[atoe], tolerance=1e-12)
  } else {
    # compare sienacpp results
    expected <- expected$sienafit[[1]]
    actual <- actual$sienafit[[1]]
    print(data.frame(targets=abs(expected$targets - actual$targets)))
    RUnit::checkEquals(expected$targets, actual$targets, tolerance=1e-12)
    print(data.frame(theta=abs(expected$theta - actual$theta)))
    RUnit::checkEquals(expected$theta, actual$theta, tolerance=1e-12)
  }
}
