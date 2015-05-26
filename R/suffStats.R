# Function to extract first and second moments of each clsTurnRes cluster
# alongside other useful stats such as % of examples in cluster and dimensions etc.
.nectr.getSuffStats <- function(x, components=NULL) {
    if(class(x) == "cTurn")
        return(.nectr.getSuffStats.cTurn(x, components))
    if(class(x) == "clsMR")
        return(.nectr.getSuffStats.clsMR(x, components))
    else
        stop("getSuffStats is defined only for cTurn and clsMR objects!")
}

.nectr.getSuffStats.cTurn <- function(x, components=NULL) {
    
    if(!is.null(components)) {
        if(length(components) > x$k) stop("Too many components specified.")
        if(!is.integer(components)) {if(!identical(components, round(components))) 
            stop("Components must be integers")}
        if(min(components) < 1 || max(components) > k) stop("Components out of range")
        k <- length(components)
    } else {
        k <- x$k
        components <- 1:k
    }
    
    dataset <- .nectr.getData(x, append.cls = FALSE, omit.na = TRUE)
    d <- ncol(dataset)
    mu <- matrix(0, k, d)
    S <- lapply(1:k, function(dummy) diag(d))
    pi <- numeric(k)
    
    cls.num <- x$cluster[!is.na(x$cluster)]
    for(i in components) {
        rws <- cls.num == i
        if(sum(rws)>1) {
            sub <- dataset[rws,]
            mu[i,] <- colMeans(sub)
            S[[i]] <- cov(sub)
            dimnames(S[[i]]) <- NULL
            pi[i] <- sum(rws)
        } else {
            warning(sum(rws), " datapoints associated with component ", i, ". Component skipped.")
            k <- k-1
            mu <- mu[-i,]
            S[i] <- NULL
            pi <- pi[-i]
        }
    }
    
    pi <- pi/sum(pi)
    n <- length(x$cluster)      # Dataset excludes non-clustered points. This vector includes all points.
    out <- list(pi=pi, mu=mu, S=S, k=k, d=d, n=n, dsn=x$dataset)
    class(out) <- ".nectr.suffStats"
    return(out)
}

.nectr.getSuffStats.clsMR <- function(x, components=NULL) {
    
    if(is.null(components)) stop("Components must be specified when using clsMR objects.")
    if(is.null(x$keep)) stop("Unable to process clsMR object. Keep = TRUE must be specified\n",
                             "in order to extract clusters for GMM.")
    
    if(length(components) > nrow(x$agglom)) stop("Too many components specified.")
    if(!is.integer(components)) {if(!identical(components, round(components))) 
        stop("Components must be integers")}
    if(min(components) < 1 || max(components) > nrow(x$agglom)) stop("Components out of range")
    
    k <- length(components)
    dataset <- .nectr.getData(x, append.cls = FALSE, omit.na = FALSE)
    d <- ncol(dataset)
    n <- nrow(dataset)
    mu <- matrix(0, k, d)
    S <- lapply(1:k, function(dummy) diag(d))
    pi <- numeric(k)
    
    
    # GET CLUSTER ASSIGNMENTS
    rowWeights <- matrix(0, n, k)
    pos <- cumsum(x$turndata[ ,"k"])         # cluster cutoffs in the clsMR hierarchy
    for(i in 1:k) {
        cl <- components[i]
        lvl <- which.max(cl <= pos)          # Level = hierarchy position (which TURN-RES run)
        local.k <- cl - c(0,pos)[lvl]        # Local k = cluster number within that level.
        rowWeights[ ,i] <- ifelse(is.na(x$keep[ ,lvl]), 0, x$keep[ ,lvl] == local.k)
    }
    rSums <- pmax(rowSums(rowWeights), 1)    # Shift all zero sum rows -> 1 for division purposes
    rowWeights <- rowWeights / matrix(rSums, n, k)
    
    # > Row Weights are to avoid a user selecting a bunch of large overlapping clusters
    # > which would skew all of the sufficient stats.
    
    # CALCULATE SUFF STATS
    for(i in 1:k) {
        rws <- rowWeights[ ,i]
        if(sum(rws)>1) {
            sub <- dataset * matrix(rws, n, d)
            mu[i, ] <- colMeans(sub)
            S[[i]] <- cov(sub)
            dimnames(S[[i]]) <- NULL
            pi[i] <- sum(rws)
        } else {
            warning(sum(rws), " datapoints associated with component ", i, ". Component skipped.")
            k <- k-1
            mu <- mu[-i,]
            S[i] <- NULL
            pi <- pi[-i]
        }
    }
    
    # TEST IF MEANS ARE TOO CLOSE
    # (impossible to gauge scale if <=2 clusters)
    if(k >2) {
        scale <- apply(mu, 2, max) - apply(mu, 2, min)
        min.dist <- sqrt(sum((scale/5)^2))
        for(uK in 1:(k-1)) {
            for(lK in (uK+1):k) {
                if(sqrt(sum((mu[uK,] - mu[lK,])^2)) < min.dist) {
                    message("Means of cluster ", components[uK], " and cluster ", components[lK],
                            " are close relative to the variation in the dataset.")
                    cat("Cluster ",components[uK]," mean = [", paste(mu[uK,], collapse=", "), "]\n")
                    cat("Cluster ",components[lK]," mean = [", paste(mu[lK,], collapse=", "), "]\n")
                    lin <- NA
                    while(is.na(lin) || !grep("^(y|yes|n|no)", lin, ignore.case=TRUE)){
                        lin <- readline("Do you wish to proceed? [y/n]")
                    }
                    if(length(grep("^(n|no)", lin, ignore.case=TRUE))>0) stop("User aborted procedure")
                }
            }
        }
    }
    
    pi <- pi/sum(pi)
    out <- list(pi=pi, mu=mu, S=S, k=k, d=d, n=n, dsn=x$dataset)
    class(out) <- ".nectr.suffStats"
    return(out)
}