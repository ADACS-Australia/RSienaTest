# Function to substitute all occurunces of the variable 'x' in an expression
substituteX <- function(x, expr) {
    do.call("substitute", list(expr, list(x = x)))
}

clusterEvalQ.SplitByRow <- function(cl, expr, xgrid) {
    # Split the xgrid array into subsets/batches
    xbatches <- snow::splitRows(xgrid, length(cl))

    # Create sets of expressions, replacing 'x' in each expression
    # with the correct subset of the xgrid array
    exprs <- lapply(xbatches, substituteX, expr = substitute(expr))

    # Evaluate the sets of expressions on each worker
    snow::docall(c, snow::clusterApply(cl, exprs, eval))
}

clusterExport.mpi.fast <- local({
    env <- as.environment(1) ## .GlobalEnv
    gets <- function(n, v) {
        assign(n, v, envir = env)
        NULL
    }
    function(cl, list, envir = .GlobalEnv) {
        ## do this with only one clusterCall--loop on workers?
        for (name in list) {
            clusterCall.mpi.fast(cl, gets, name, get(name, envir = envir))
        }
    }
})

clusterCall.mpi.fast <- function(cl, fun, ...) {
    snow::checkCluster(cl)

    # create packet and serialise
    value <- list(fun = fun, args = list(...), return = TRUE, tag = NULL)
    data_ <- list(type = "EXEC", data = value, tag = NULL)
    data <- serialize(data_, NULL)

    # send packet to all workers
    for (i in seq(along = cl)) {
        node <- cl[[i]]
        Rmpi::mpi.isend(x = data, type = 4, dest = node$rank, tag = node$SENDTAG, comm = node$comm)
    }

    # recieve results from workers and return
    snow::checkForRemoteErrors(lapply(cl, snow::recvResult))
}
