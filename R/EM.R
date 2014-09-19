   # dyn.load("C:/Users/pc/Documents/pkg_dev/nectr/src/rowWhichMaxC.dll")
getEMGPs <- function(data, centers, cls.prob = NULL, move.restrict = NA, 
                     eps.target = 0.1, trunc = 0.5, max.iter = 100, silent = FALSE) {
    
    #Set up
    #-----------
    require(mvtnorm)
    m = nrow(data)
    n = ncol(data)
    mu <- centers
    if(!(is.list(mu) & !is.data.frame(mu))) stop("Centers must be a list of vectors.")
    if(!identical(unname(mu), unique(mu))) stop("Centers have repeated values. They must be distinct")
    k = length(mu)
    S <- list(); for(i in 1:k) S[[i]] <- diag(n)
    w <- matrix(0, m, k)
    phi <- vector()
    if(is.null(cls.prob)) {
        for(i in 1:k) phi[i] <- 1/k
    } else {
        if(is.numeric(cls.prob) & length(cls.prob) == k) {
            phi <- cls.prob/sum(cls.prob)
        } else {
            stop("cls.prob must be vector with same dimension as centers")
        }
    }
    if(!is.na(move.restrict) & !(is.numeric(move.restrict) & length(move.restrict)==1))
        stop("move.restrict must be a scalar value")
    llh.history <- vector(mode="numeric", length=max.iter+1)
    llh.history[1] <- -Inf
    if(is.data.frame(data)) data <- as.matrix(data)
    
    #EM Loop
    for(cnt in 1:max.iter) {
        #E: Set-up
        mu <- lapply(mu, as.matrix)
        
        #E: Assign un-normalised wj's
        for(j in 1:k) {
            w[ ,j] <- dmvnorm(data, mu[[j]], S[[j]])*phi[j]
        }
        #E: Normalise wj's
        w.maxs <- .Call("rowWhichMaxC", w)       #[,1] = rowMax, [,2] = which.max
        zeros <- w.maxs[, 1] < trunc/m
        w <- w/rowSums(w)
        for(i in 1:k) w[zeros,i] <- 0
        
        #Convergence
        llh <- 0
        for(j in 1:k) {
            llh <- llh + sum(log(dmvnorm(data[w.maxs[ ,2] == j, ], mu[[j]], S[[j]])))
        }
        eps <- llh - llh.history[cnt]
        if(!silent) cat("EM Iteration ", cnt,": ", as.numeric(eps),"\n")
        llh.history[cnt+1] <- llh
        if(abs(eps) < eps.target) break
        
        #M-step
        #--------------
        
        sum.w <- colSums(w)
        phi <- sum.w/m
        
        for(j in 1:k) {
            data.centered <- sqrt(w[,j]) * (data - matrix(mu[[j]], m, n, byrow=T))
            sigma.sum <- t(data.centered) %*% data.centered
            S[[j]] <- sigma.sum/sum.w[j]
        }
        
        #Update mu, but impose restrictions (if any)
        for(i in 1:k) {
            move <- colSums(w[ ,i] * data)/sum.w[i] - mu[[i]]
            if(is.na(move.restrict)) {
                mu[[i]] <- mu[[i]] + move
            } else {
                len <- sqrt(sum(move^2))
                mu[[i]] <- mu[[i]] + move*pmin(move.restrict,len)/len
            }
        }
    }
    
    if(abs(eps) >= eps.target) warning("EM Algorithm did not converge within max.iter")
    for(i in 1:k) dimnames(mu[[i]]) <- NULL
    for(i in 1:k) dimnames(S[[i]]) <- NULL
    out <- list()
    out$llh.history <- llh.history[2:(cnt+1)]
    out$phi <- phi
    out$mu <- mu
    out$sigma <- S
    return(out)
}



getEMClusters <- function(data, mu = NA, sigma = NA, phi = NA, params = NA) {
    
    #Set up
    #-----------
    require(mvtnorm)
    m = nrow(data)
    n = ncol(data)
    if(is.list(params) &!is.data.frame(params)) {
        if(!(is.list(params) & !is.data.frame(params))) stop("params supplied must be a list.")
        mu <- tryCatch(params[["mu"]], error = function(e) e)
        if(inherits(mu, "error")) stop("Params list does not contain 'mu' slot.")
        sigma <- tryCatch(params[["sigma"]], error = function(e) e)
        if(inherits(sigma, "error")) stop("Params list does not contain 'sigma' slot.")
        phi <- tryCatch(params[["phi"]], error = function(e) e)
        if(inherits(phi, "error")) stop("Params list does not contain 'phi' slot.")
    } else {
        if(!is.na(params)) warning("Params is not a list and will be ignored")
        if(is.na(mu)) stop("mu is missing. Either supply 'params' or 'mu', 'sigma' and 'phi'.")
        if(is.na(sigma)) stop("sigma is missing. Either supply 'params' or 'mu', 'sigma' and 'phi'.")
        if(is.na(phi)) stop("phi is missing. Either supply 'params' or 'mu', 'sigma' and 'phi'.")
    }
    if(!(is.list(mu) & !is.data.frame(mu))) stop("Centers must be a list of vectors.")
    if(!identical(mu, unique(mu))) stop("Some mu(j) are the same. They must be distinct.")
    k = length(mu)
    if(!(is.list(sigma) & !is.data.frame(sigma))) stop("Sigma must be a list of matrices.")
    if(length(sigma) != k) stop("Different number of covariance matrices to means.")
    for(i in 1:k) {
        if(is.data.frame(sigma[[i]])) sigma[[i]] <- as.matrix(sigma[[i]])
        if(!is.matrix(sigma[[i]])) stop("Sigma must be a list of matrices")
        if(nrow(sigma[[i]]) != ncol(sigma[[i]])) stop("Covariance mx ", i, " not square.")
        if(nrow(sigma[[i]]) != n) stop("Covariance mx ", i, " of different dimension to mu.")
    }
    
    w <- matrix(0, m, k)
    if(is.data.frame(data)) data <- as.matrix(data)
    if(!is.matrix(data)) stop("data must be matrix or data.frame.")
    if(!is.numeric(data)) stop("data must be numeric only.")
    mu <- lapply(mu, as.matrix)
    
    #Assign un-normalised wj's
    for(j in 1:k) {
        w[ ,j] <- dmvnorm(data, mu[[j]], sigma[[j]])*phi[j]
    }
    #Normalise wj's
    #w <- w/rowSums(w)
    
    #Which.max
    w <- .Call("rowWhichMaxC", w)       #[,1] = rowMax, [,2] = which.max
    w <- data.frame(cluster = w[ ,2], prob = w[ ,1])
    return(w)
}

formaliseClusters <- function(x, clusters = NULL, data = NULL, ...) {
    
    Call <- match.call(expand.dots = TRUE)
    if(class(x) == "cTurn") {
        if(is.null(clusters)) clusters <- as.numeric(x$table$res)
        centr <- as.list(as.data.frame(t(x$means[clusters, ])))
    } else if(class(x) == "clsMR") {
        if(is.null(clusters)) stop("if passing clsMR object, clusters must be specified")
        centr <- as.list(as.data.frame(t(x$means[clusters, ])))
    } else {
        stop("object class '", class(x),"' not recognised. Supply either cTurn or clsMR object.")
    }        
            
    attr(centr, "names") <- NULL
    if(is.null(data)) {
        e <- try(data <- get(x$dataset, envir = .GlobalEnv), silent = TRUE)
        if(inherits(e, "try-error")) {
            stop("Error obtaining data from $dataset.name - unable to find in workspace.\n",
                 "Original dataset can be passed in using data = ... argument")
        }
    }
    
    if(is.null(Call[["max.iter"]])) max.iter <- 100 else max.iter <- Call[["max.iter"]]
    if(is.null(Call[["eps.target"]])) eps.target <- 0.1 else max.iter <- Call[["eps.target"]]
    message("Performing EM ascent on dataset - this may take some time.\nCurrent params are: ",
            "max iterations:", max.iter, ", convergence target: ", eps.target)
    params <- getEMGPs(data, centr, ...)
    message("Parameters estimated. Scoring all data...")
    cls <- getEMClusters(data, params = params)
    params$cluster <- cls$cluster
    params$cluster.prob <- cls$prob
    return(params)
}

