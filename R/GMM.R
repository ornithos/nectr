# dyn.load("C:/Users/pc/Documents/pkg_dev/nectr/src/rowWhichMaxC.dll")



.nectr.fitGMM <- function(params, hyperparams, eps.stop = 0.01, max.iter = 100, silent = FALSE) {
    
    #Set up
    #-----------
    .nectr.checkParams(params)
    data <- .nectr.getData(params$dsn)
    if(is.data.frame(data)) data <- as.matrix(data)
    
    # ..........................................................................
    # Calculate the density of Gaussian with mean and cov of each datapoint
    # in x. Credit belongs entirely to mvtnorm package. Since it's for internal
    # use, there is no need of the error catching in dmvnorm.
    # 
    # ** Shamelessly stolen from code by Friedrich Leisch and Fabian Scheipl **
    .internal.dmvnorm <- function(x, mean, cov, log_units = FALSE) {
        dec <- tryCatch(chol(cov), error = function(e) e)
        if (inherits(dec, "error")) {
            warning("Covariance matrix is not positive definite")
            browser()
            logdensity <- rep.int(-Inf, x)
        }
        else {
            tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
            rss <- colSums(tmp^2)
            logdensity <- -sum(log(diag(dec))) - 0.5 * ncol(x) * log(2 * pi) - 0.5 * rss
        }
        if(log_units) return(logdensity)
        return(exp(logdensity))
    }
    
    # ..........................................................................
    
    # Trace - to make llh more readable
    Tr <- function(x) sum(diag(x))
    
    # Function to Calc (Log) Likelihood
    calc.llh <- function(par, hyper, pri, pi, mu, S) {
        llh <- numeric(4)
        
        # DATA TERMS: due to the sum inside log, all densities must be precalcd to take
        # advantage of linear algebra routines. Waste of memory, but needs -> C to do better.
        data_term_store <- matrix(0, par$n, par$k)
        for(cK in 1:par$k) {
            data_term_store[,cK] <- .internal.dmvnorm(data, mu[cK,], S[[cK]], log_units=FALSE) * pi[cK]
        }
        llh[1] <- sum(log(rowSums(data_term_store)))

        # PRIOR TERMS: these are a simple sum and so can be looped over in the usual way
        for(cK in 1:par$k) {
            llh[2] <- llh[2] + (hyper$alpha[cK] - 1)*log(pi[cK])
            llh[3] <- llh[3] + .internal.dmvnorm(matrix(mu[cK,],ncol=par$d), pri$mu[cK,], S[[cK]] /
                                               hyper$beta[cK], log_units=TRUE)
            llh[4] <- llh[4] - 0.5*(hyper$nu[cK]+par$d+1)*log(det(S[[cK]]))
            llh[4] <- llh[4] - 0.5*hyper$nu[cK]*Tr(pri$S[[cK]] %*% solve(S[[cK]]))
        }
        return(c(sum(llh),llh))
    }
    # ..........................................................................
    
    
    params$r <- with(params, matrix(0, n, k))
    if(is.data.frame(params$mu)) within(params, mu <-  as.matrix(mu))
    
    gauss.ldensity <- with(params, matrix(0, n, k))
    Nk <- numeric(params$k)
    prior <- params
    llh.history <- matrix(0,max.iter+1,5)
    llh.history[1,] <- rep(-Inf,5)
    eps <- 0
    
	# data = n * d
	# number of clusters = k;  score of each datapoint in each cluster = r.
	# 
    #EM Loop
    out <- with(params, {
        for(cI in 1:max.iter) {    
        #E: Assign r_nk's
        #-----------------
            for(cK in 1:k) {
                gauss.ldensity[ ,cK] <- .internal.dmvnorm(data, mu[cK,], S[[cK]], log_units=TRUE)
                r[ ,cK] <- gauss.ldensity[ ,cK] + log(pi[cK])
            }
            r <- exp(r)/rowSums(exp(r))
            Nk <- colSums(r)

            # Test for Convergence
            llh <- calc.llh(params, hyperparams, prior, pi, mu, S)
            eps <- llh[1] - llh.history[cI,1]
            if(!silent) cat("EM Iteration ", cI,": log likelihood change: ", as.numeric(eps),"\n")
            llh.history[cI+1,] <- llh
            if(abs(eps) < eps.stop) break
        
        #M-step
        #--------------
            # MIXING PROBABILITIES, PI
            if(all(hyperparams$infinities$alpha == FALSE)) {
                for(cK in 1:k) pi[cK] <- hyperparams$alpha[cK] + Nk[cK] - 1
                pi <- pi / sum(pi)
            }
        
            # COMPONENT MEANS, MU
            if(all(hyperparams$infinities$beta == FALSE)) {
                for(cK in 1:k) mu[cK, ] <- colSums(r[ ,cK]*data) + hyperparams$beta[cK]*prior$mu[cK, ]
                mu <- mu / matrix(Nk + hyperparams$beta, k, d)
            }
            
            # COMPONENT COVARIANCE, SIGMA
            if(all(hyperparams$infinities$nu == FALSE & hyperparams$infinities$beta == FALSE)) {
                for(cK in 1:k) {
                    cov.tmp <- t(Reduce("-", list(t(data), mu[cK,])))
                    prior.outer <- outer(mu[cK,] - prior$mu[cK,], mu[cK,] - prior$mu[cK,])
                    S[[cK]] <- t(r[ ,cK]*cov.tmp) %*% cov.tmp
                    S[[cK]] <- S[[cK]] + hyperparams$beta[cK]*prior.outer + hyperparams$nu[cK]*prior$S[[cK]]
                    S[[cK]] <- S[[cK]] / (Nk[cK] + hyperparams$nu[cK] + d + 2)
                }
            }
            
        }
        retval <- list(pi=pi, mu=mu, S=S, k=k, d=d, n=n, dsn=dsn)
        retval$cluster <- .Call("rowWhichMaxC", r)[ ,2]    #[,1] = rowMax, [,2] = which.max
        retval$llh_ascent <- llh.history[2:(cI+1), 1]
        retval$monotone <- !is.unsorted(retval$llh_ascent)
        retval$iter <- cI
        retval
    })
    if(abs(out$llh_ascent[out$iter] - out$llh_ascent[out$iter-1]) >= eps.stop) 
        warning("EM Algorithm did not converge within max.iter")
    if(!out$monotone) warning("Log Likelihood ascension was not monotone.")
    class(out) <- ".nectr.GMM"
    return(out)
}



scoreGMM <- function(data, mu = NA, sigma = NA, phi = NA, params = NA) {
    
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
    if(is.null(Call[["eps.stop"]])) eps.stop <- 0.1 else max.iter <- Call[["eps.stop"]]
    message("Performing EM ascent on dataset - this may take some time.\nCurrent params are: ",
            "max iterations:", max.iter, ", convergence target: ", eps.stop)
    params <- fitGMM(data, centr, ...)
    message("Parameters estimated. Scoring all data...")
    cls <- scoreGMM(data, params = params)
    params$cluster <- cls$cluster
    params$cluster.prob <- cls$prob
    return(params)
}

