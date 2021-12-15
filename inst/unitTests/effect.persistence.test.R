# This file contains the effect before/after style tests. Have a look at
# R/runit.R for how to run them.
source('asserts.R')

# currently untested effect groups
#
#   behaviorBipBipObjective
#   behaviorBipartiteObjective
#   behaviorOneModeSymObjective
#   behaviorOneOneModeObjective
#   behaviorSymSymObjective
#   behaviorSymmetricObjective
#   bipartiteBipartiteObjective
#   bipartiteNonSymmetricObjective
#   bipartiteObjective
#   bipartiteSymmetricObjective
#   covarABNetNetObjective
#   covarABehaviorBipartiteObjective
#   covarANetNetObjective
#   covarBBehaviorBipartiteObjective
#   covarBNetNetObjective
#   covarBehaviorNetObjective
#   covarBehaviorObjective
#   covarBipartiteObjective
#   dyadBipartiteObjective
#   dyadObjective
#   nonSymmetricBipartiteObjective
#   nonSymmetricNonSymmetricObjective
#
#   unspecifiedBehaviorInteraction
#   unspecifiedNetInteraction
#
#   behaviorBipartiteRate
#   behaviorOneModeRate
#   behaviorRate
#   behaviorSymmetricRate
#   bipartiteRate
#   covarBehaviorOneModeRate
#   covarBehaviorRate
#   covarBipartiteRate
#   covarNonSymmetricRate
#   covarSymmetricRate
#   nonSymmetricRate
#   symmetricRate

#' estimates a model with fixed seed and fast algorithm settings
#' @param data siena data
#' @param eff siena effects
run_effect_test_model <- function(data, eff) {
  model <- sienaModelCreate(projname='effect_test', seed=12345, nsub=1, n3=3)
  siena07(model, data=data, effects=eff, batch=T)
}

#' @param group effect group name
#' @return vector with the shortNames of all evaluation effects in group
short_names <- function(effect_info_filter) {
  select <- apply(RSienaTest::allEffects, 1, effect_info_filter)
  sort(levels(factor(RSienaTest::allEffects[which(select),]$shortName)))
}

#' @param data siena data
#' @param group effect group name
#' @param ... passed along to includeEffects
#' @return siena effect object with all effects in the group
all_non_gmm_effects <- function(data, group, ...) {
  eff <- getEffects(data)
  filter <- function(effect_info) {
    effect_info['type'] != 'gmm' && effect_info['effectGroup'] == group
  }
  for (name in short_names(filter)) {
    eff <- includeEffects(eff, name, ..., character=T)
  }
  eff
}

test.effect.symmetric.objective <- function() {
  group <- 'symmetricObjective'
  if (skip_recording(group)) return()
  dn <- textConnection(NULL, 'w')
  sink(dn)
  net <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
  data <- sienaDataCreate(net)
  eff <- all_non_gmm_effects(data, group)
  ans <- run_effect_test_model(data, eff)
  sink()
  close(dn)
  check_model_persistence(group, ans)
}

test.effect.non.symmetric.objective <- function() {
  group <- 'nonSymmetricObjective'
  if (skip_recording(group)) return()
  dn <- textConnection(NULL, 'w')
  sink(dn)
  net <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
  data <- sienaDataCreate(net)
  eff <- all_non_gmm_effects(data, group)
  ans <- run_effect_test_model(data, eff)
  sink()
  close(dn)
  check_model_persistence(group, ans)
}

test.effect.covar.symmetric.objective <- function() {
  group <- 'covarSymmetricObjective'
  if (skip_recording(group)) return()
  dn <- textConnection(NULL, 'w')
  sink(dn)
  net <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
  atr <- varCovar(s50a)
  data <- sienaDataCreate(net, atr)
  eff <- all_non_gmm_effects(data, group, interaction1='atr')
  ans <- run_effect_test_model(data, eff)
  sink()
  close(dn)
  check_model_persistence(group, ans)
}

test.effect.covar.non.symmetric.objective <- function() {
  group <- 'covarNonSymmetricObjective'
  if (skip_recording(group)) return()
  dn <- textConnection(NULL, 'w')
  sink(dn)
  net <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
  atr <- varCovar(s50a)
  data <- sienaDataCreate(net, atr)
  eff <- all_non_gmm_effects(data, group, interaction1='atr')
  ans <- run_effect_test_model(data, eff)
  sink()
  close(dn)
  check_model_persistence(group, ans)
}

test.effect.behavior.objective <- function() {
  group <- 'behaviorObjective'
  if (skip_recording(group)) return()
  dn <- textConnection(NULL, 'w')
  sink(dn)
  atr <- sienaDependent(s50a, type='behavior')
  data <- sienaDataCreate(atr)
  eff <- all_non_gmm_effects(data, group)
  ans <- run_effect_test_model(data, eff)
  sink()
  close(dn)
  check_model_persistence(group, ans)
}

ignored.effect.behavior.one.mode.objective <- function() {
  group <- 'behaviorOneModeObjective'
  if (skip_recording(group)) return()
  dn <- textConnection(NULL, 'w')
  sink(dn)
  net <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
  atr <- sienaDependent(s50a, type='behavior')
  data <- sienaDataCreate(net, atr)
  eff <- all_non_gmm_effects(data, group, interaction1='net')
  ans <- run_effect_test_model(data, eff)
  sink()
  close(dn)
  check_model_persistence(group, ans)
}

ignored.effect.objective.gmm.onemode <- function() {
  # sink(NULL)
  group <- 'gmm.nonSymmetricObjective'
  if (skip_recording(group)) return()
  net <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
  data <- sienaDataCreate(net)
  # add gmm effects (automatic is not a good idea, need to be align to the objective)
  eff <- getEffects(data)
  eff <- includeEffects(eff, density)
  eff <- includeEffects(eff, density, type='gmm')
  eff <- includeEffects(eff, recip)
  eff <- includeEffects(eff, newrecip, type='gmm')
  eff <- includeEffects(eff, realrecip, type='gmm')
  eff <- includeEffects(eff, persistrecip, type='gmm')
  eff <- includeEffects(eff, transTrip)
  eff <- includeEffects(eff, transTrip, type='gmm')
  eff <- includeEffects(eff, realtrans, type='gmm')
  # run gmm model
  model <- sienaModelCreate(seed=12345, nsub=1, n3=13, dolby=F)
  ans <- sienacpp(model, data=data, effects=eff, logLevelConsole='INFO')
  # sink()
  check_model_persistence(group, ans)
}

# tried the next 3 but I guess with to simple model/data, estimation finishes with errors
ignored.effect.non.symmetric.symmetric.objective <- function() {
  group <- 'nonSymmetricSymmetricObjective'
  if (skip_recording(group)) return()
  # sink(NULL)
  net1 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
  net2 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
  data <- sienaDataCreate(net1, net2)
  eff <- all_non_gmm_effects(data, group, interaction1='net2')
  ans <- run_effect_test_model(data, eff)
  # sink()
  check_model_persistence(group, ans)
}

ignored.effect.covar.net.net.objective <- function() {
  group <- 'covarNetNetObjective'
  if (skip_recording(group)) return()
  dn <- textConnection(NULL, 'w')
  sink(dn)
  net1 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
  net2 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
  atr <- varCovar(s50a)
  data <- sienaDataCreate(net1, net2, atr)
  eff <- all_non_gmm_effects(data, group, interaction1='atr', interaction2='net2')
  ans <- run_effect_test_model(data, eff)
  sink()
  close(dn)
  check_model_persistence(group, ans)
}

ignored.effect.triple.network.objective <- function() {
  group <- 'tripleNetworkObjective'
  if (skip_recording(group)) return()
  dn <- textConnection(NULL, 'w')
  sink(dn)
  net1 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
  net2 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
  net3 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
  data <- sienaDataCreate(net1, net2, net3)
  eff <- all_non_gmm_effects(data, group, interaction1='net2', interaction2='net3')
  ans <- run_effect_test_model(data, eff)
  sink()
  close(dn)
  check_model_persistence(group, ans)
}
