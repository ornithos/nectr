.nectr.getData <- function(x, append.cls = FALSE, omit.na = FALSE) {
    dat <- get(x$dataset, envir = .GlobalEnv)
    if(omit.na) {
        rws <- !is.na(x$cluster)
    } else {
        rws <- 1:nrow(dat)
    }
    
    if(append.cls) {
        return(cbind(dat[rws], cls = x$cluster[rws]))
    } else {
        return(dat[rws,])
    }
}

# Function to extract first and second moments of each clsTurnRes cluster
# alongside other useful stats such as % of examples in cluster and dimensions etc.
.nectr.getSuffStats <- function(x) {
    if(class(x) != "cTurn") stop("getSuffStats is defined only for cTurn objects!")
    k <- x$k
    dataset <- .nectr.getData(x, append.cls = FALSE, omit.na = TRUE)
    d <- ncol(dataset)
    n <- nrow(dataset)
    mu <- matrix(0, k, d)
    S <- lapply(1:k, function(dummy) diag(d))
    pi <- numeric(k)
    
    cls.num <- x$cluster[!is.na(x$cluster)]
    for(i in 1:k) {
       rws <- cls.num == i
       sub <- dataset[rws,]
       mu[i,] <- colMeans(sub)
       S[[i]] <- cov(sub)
       pi[i] <- sum(rws)/n
    }
    
    out <- list(pi=pi, mu=mu, S=S, k=k, d=d)
    class(out) <- ".cTurn.suffStats"
    return(out)
}

.nectr.checkParams <- function(x) {
    err.msg <- "Object is not of valid GMM parameter type"
    if(class(x) != ".cTurn.suffStats") stop(err.msg, ": (class error)")
    slots <- names(x)
    if(!(length(slots) == 5 & length(union(slots, c("pi","mu","S","k","d"))==5))) {
        stop(err.msg, ": (incorrect object structure)")}
    if(sum(x$pi) > 1 + 10^-9 | sum(x$pi) < 1 - 10^-9) stop(err.msg, ": (pi doesn't sum to 1)")
    if(!(length(x$k) == 1 & length(x$d) == 1 & is.numeric(x$k) & is.numeric(x$d))) {
        stop(err.msg, ": (missing num clusters k or dimension d)")}
    if(!all(dim(x$mu) == c(x$k, x$d))) {
        stop(err.msg, ": (mu table is wrong dimensions - expecting ", x$k, " by ", x$d, ".)")}
    if(!(is.list(x$S) & !is.data.frame(x$S) & length(x$S) == x$k)) {
        stop(err.msg, ": (S is not a list of ", x$k, " elements)")}

    for(i in 1:x$k) {
        if(!is.matrix(x$S[[i]])) stop(err.msg, ": (S[", i, "] is not a matrix)")
        if(!all(dim(x$S[[i]]) == c(x$d,x$d))) stop(err.msg, ": (S[", i, "] is not d x d)")
    }
    return(invisible(TRUE))
}