# This file contains the model options before/after style tests. Have a look
# at R/runit.R for how to run them.
source('asserts.R')

# tests generator
test_environment <- environment()
model_option <- function(name, values, ...) {
  for (value in values) {
    (function() {
      op <- list(projname='model_test', seed=12345, nsub=3, n3=500, ...)
      op[[name]] <- value
      model <- do.call('sienaModelCreate', op)
      fn_name <- paste('test.model', name, value, paste(names(list(...)), sep='_'), sep='.')
      assign(fn_name, function() {
        if (skip_recording(fn_name)) return()
        dn <- textConnection(NULL, 'w')
        sink(dn)
        net <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
        data <- sienaDataCreate(net)
        eff <- getEffects(data)
        eff <- includeEffects(eff, transTrip)
        ans <- siena07(model, data=data, effects=eff, batch=T)
        sink()
        close(dn)
        check_model_persistence(fn_name, ans)
      }, envir=test_environment)
    })()
  }
}

model_option('dolby', c(T,F))
model_option('maxlike', T)
model_option('cond', T)
model_option('findiff', T)
model_option('findiff', T, cond=T)
