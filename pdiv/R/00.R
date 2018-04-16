
## apply function that uses parallelization if available

.parallel <-
    (Sys.info()["sysname"] != "Windows") &&
    requireNamespace("parallel")

.mysapply <-
    if(.parallel) {
        function(x, fun)
            unlist(parallel::mclapply(x, fun, mc.cores = parallel::detectCores()))
    } else
        function(x, fun)
            sapply(x, fun)

.mylapply <-
    if(.parallel) {
        function(x, fun)
            parallel::mclapply(x, fun, mc.cores = parallel::detectCores())
    } else
        function(x, fun)
            lapply(x, fun)

## define population covariance
.covp <- function(x, y) {
    sum((x - mean(x)) * (y - mean(y))) / length(x);
}

## convert all factors to numeric values
.columnsToNum <- function(d)
{
    for(v in names(d))
        d[[v]] <- as.numeric(as.factor( d[[v]] ))
    d
}
