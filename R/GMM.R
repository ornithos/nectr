# dyn.load("C:/Users/pc/Documents/pkg_dev/nectr/src/rowWhichMaxC.dll")

# ..........................................................................
# Calculate the density of Gaussian with mean and cov of each datapoint
# in x. Credit belongs entirely to mvtnorm package. Since it's for internal
# use, there is no need of the error catching in dmvnorm.
# 
# ** Shamelessly stolen from code by Friedrich Leisch and Fabian Scheipl **
.nectr.dmvnorm <- function(x, mean, cov, log_units = FALSE) {
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


.nectr.fitGMM <- function(params, hyperparams, eps.stop = 0.001, max.iter = 100, 
                          skip.llh = 1, silent = FALSE) {
    
    #Set up
    #-----------
    .nectr.checkParams(params)
    .nectr.checkParams(hyperparams)
    data <- .nectr.getData(params$dsn)
    if(is.data.frame(data)) data <- as.matrix(data)
      
    # Trace - to make llh more readable
    Tr <- function(x) sum(diag(x))
    
    # Function to Calc (Log) Likelihood
    calc.llh <- function(par, hyper, pri, pi, mu, S) {
        llh <- numeric(4)
        # DATA TERMS: due to the sum inside log, all densities must be precalcd to take
        # advantage of linear algebra routines. Waste of memory, but needs -> C to do better.
        data_term_store <- matrix(0, par$n, par$k)
        for(cK in 1:par$k) {
            data_term_store[,cK] <- .nectr.dmvnorm(data, mu[cK,], S[[cK]], log_units=FALSE) * pi[cK]
        }
        llh[1] <- sum(log(rowSums(data_term_store)))

        # PRIOR TERMS: these are a simple sum and so can be looped over in the usual way
        for(cK in 1:par$k) {
            if(is.finite(hyper$alpha[cK])) llh[2] <- llh[2] + (hyper$alpha[cK] - 1)*log(pi[cK])
            if(is.finite(hyper$beta[cK])) {
                llh[3] <- llh[3] + .nectr.dmvnorm(matrix(mu[cK,],ncol=par$d), pri$mu[cK,], S[[cK]] /
                                                   hyper$beta[cK], log_units=TRUE)
            }
            if(is.finite(hyper$nu[cK])) {
                llh[4] <- llh[4] - 0.5*(hyper$nu[cK]+par$d+1)*log(det(S[[cK]]))
                llh[4] <- llh[4] - 0.5*hyper$nu[cK]*Tr(pri$S[[cK]] %*% solve(S[[cK]]))
            }
        }
        return(c(sum(llh),llh))
    }
    # ..........................................................................
    
    
    if(is.data.frame(params$mu)) within(params, mu <-  as.matrix(mu))
    
    llh.history <- matrix(0,max.iter+1,5)
    llh.history[1,] <- rep(-Inf,5)
    eps <- 0
    
	# data = n * d
	# number of clusters = k;  score of each datapoint in each cluster = r.
	# 

    k <- params$k
    n <- params$n
    d <- params$d
    p <- params$pi
    mu <- params$mu
    S <- params$S
    prior <- params
    
    gauss.ldensity <- matrix(0, n, k)
    Nk <- numeric(k)
    r <- matrix(0, n, k)
    
    #EM Loop
    
    for(cI in 1:max.iter) {    
    #E: Assign r_nk's
    #-----------------
        for(cK in 1:k) {
            r[ ,cK] <- .nectr.dmvnorm(data, mu[cK,], S[[cK]])*p[cK]
        }
        r <- r/rowSums(r)
        Nk <- colSums(r)

        # Test for Convergence
        if(cI %% skip.llh == 0) {
            llh <- calc.llh(params, hyperparams, prior, p, mu, S)
            eps <- (llh[1] - llh.history[(cI %/% skip.llh),1])/skip.llh
            if(!silent) cat("\r","EM Iteration ", cI,": log likelihood change: ", as.numeric(eps),"      ")
            llh.history[(cI %/% skip.llh)+1,] <- llh
            if(abs(eps) < eps.stop) break
        }
    
    #M-step
    #--------------
        # MIXING PROBABILITIES, PI
        if(all(hyperparams$infinities$alpha == FALSE)) {
            for(cK in 1:k) p[cK] <- hyperparams$alpha[cK] + Nk[cK] - 1
            p <- p / sum(p)
        }
        
        # COMPONENT MEANS, MU
        if(all(hyperparams$infinities$beta == FALSE)) {
            for(cK in 1:k) mu[cK, ] <- colSums(r[ ,cK]*data) + hyperparams$beta[cK]*prior$mu[cK, ]
            mu <- mu / matrix(Nk + hyperparams$beta, k, d)
        }
            
        # COMPONENT COVARIANCE, SIGMA
        if(all(hyperparams$infinities$nu == FALSE & hyperparams$infinities$beta == FALSE)) {
            for(cK in 1:k) {
                #cov.tmp <- t(Reduce("-", list(t(data), mu[cK,]))) # - slightly more efficient, but less
                # readable, and I believe the 2 tranpose operations negate most of the efficiency.
                cov.tmp <- data - matrix(mu[cK,], n, d, byrow=TRUE)
                prior.outer <- outer(mu[cK,] - prior$mu[cK,], mu[cK,] - prior$mu[cK,])
                S[[cK]] <- t(r[ ,cK]*cov.tmp) %*% cov.tmp
                S[[cK]] <- S[[cK]] + hyperparams$beta[cK]*prior.outer + hyperparams$nu[cK]*prior$S[[cK]]
                S[[cK]] <- S[[cK]] / (Nk[cK] + hyperparams$nu[cK] + d + 2)
            }
        }
            
    }
    
    out <- list(pi=p, mu=mu, S=S, k=k, d=d, n=n, dsn=params$dsn)
    cls <- .Call("rowWhichMaxC", r)    #[,1] = rowMax, [,2] = which.max
    out$cluster <- data.frame(cluster = cls[ ,2], prob = cls[ ,1])
    out$noise <- NA
    out$llh_ascent <- llh.history[2:((cI %/% skip.llh)+1), 1]
    out$monotone <- !is.unsorted(out$llh_ascent)
    out$iter <- c(cI, cI %/% skip.llh)
    
    if(abs(out$llh_ascent[out$iter[2]] - out$llh_ascent[out$iter[2]-1]) >= eps.stop*skip.llh) 
        warning("EM Algorithm did not converge within max.iter")
    if(!out$monotone) warning("Log Likelihood ascension was not monotone.")
    out$fit <- TRUE
    class(out) <- ".nectr.GMM"
    return(out)
}




predict.nectr.GMM <- function(object, newdata)

scoreGMM <- function(data, mu = NA, sigma = NA, phi = NA, params = NA) {
    
    #Set up
    #-----------
    m = nrow(data)
    n = ncol(data)

    
    w <- matrix(0, m, k)
    if(is.data.frame(data)) data <- as.matrix(data)
    if(!is.matrix(data)) stop("data must be matrix or data.frame.")
    if(!is.numeric(data)) stop("data must be numeric only.")
    mu <- lapply(mu, as.matrix)
    
    #Assign un-normalised wj's
    for(j in 1:k) {
        w[ ,j] <- .nectr.dmvnorm(data, mu[[j]], sigma[[j]])*pi[j]
    }
    #Normalise wj's
    #w <- w/rowSums(w)
    
    #Which.max
    w <- .Call("rowWhichMaxC", w)       #[,1] = rowMax, [,2] = which.max
    w <- data.frame(cluster = w[ ,2], prob = w[ ,1])
    return(w)
}

