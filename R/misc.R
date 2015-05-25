.nectr.getData <- function(x, append.cls = FALSE, omit.na = FALSE) {
    if(class(x) == "cTurn") {
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
    } else if(class(x) == "character") {
        if(append.cls | omit.na) {
            message(".nectr.getData: append.cls and omit.na options are unused when x is a character")}
        dat <- try(get(x, envir = .GlobalEnv), silent=TRUE)
        if(inherits(dat, "try-error")) stop("Dataset '", x,"' does not exist in .GlobalEnv")
        return(dat)
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
       dimnames(S[[i]]) <- NULL
       pi[i] <- sum(rws)/n
    }
    
    n <- length(x$cluster)
    out <- list(pi=pi, mu=mu, S=S, k=k, d=d, n=n, dsn=x$dataset)
    class(out) <- ".cTurn.suffStats"
    return(out)
}

.nectr.checkParams <- function(x) {
    err.msg <- "Object is not of valid GMM parameter type"
    if(class(x) != ".cTurn.suffStats") stop(err.msg, ": (incorrect class)")
    slots <- names(x)
    if(!(length(slots) == 7 & length(union(slots, c("pi","mu","S","k","d","n","dsn"))==7))) {
        stop(err.msg, ": (incorrect object structure)")}
    if(sum(x$pi) > 1 + 10^-9 | sum(x$pi) < 1 - 10^-9) stop(err.msg, ": (pi doesn't sum to 1)")
    if(!(length(x$k) == 1 & length(x$d) == 1 & is.numeric(x$k) & is.numeric(x$d))) {
        stop(err.msg, ": (missing num clusters k or dimension d)")}
    if(!all(dim(x$mu) == c(x$k, x$d))) {
        stop(err.msg, ": (mu table is wrong dimensions - expecting ", x$k, " by ", x$d, ".)")}
    if(!(is.list(x$S) & !is.data.frame(x$S) & length(x$S) == x$k)) {
        stop(err.msg, ": (S is not a list of ", x$k, " elements)")}
    if(!(length(x$n) == 1 & is.numeric(x$n))) stop(err.msg, ": (n must be a numeric literal)")
    
    for(i in 1:x$k) {
        if(!is.matrix(x$S[[i]])) stop(err.msg, ": (S[", i, "] is not a matrix)")
        if(!all(dim(x$S[[i]]) == c(x$d,x$d))) stop(err.msg, ": (S[", i, "] is not d x d)")
    }
    return(invisible(TRUE))
}

# User-friendly function to generate the hyperparameters. Can be expressed as only 2 literals,
# or with alpha, beta, nu specified individually for each k. Handles both (%) and (abs#)
.nectr.getHyper <- function(x, ...) {
    dots <- list(...)
    if(class(x) != ".cTurn.suffStats") stop("x must be of type .cTurn.suffStats")
    .nectr.checkParams(x)
    
    # Each parameter can be specified as a single number, or vector of k values.
    # also, can either specify the absolute number of the multiple of n to choose.
    # eg. 'alpha' is the multiplier, 'alpha_abs' is the absolute number
    processParam <- function(f_inp, param) {
        if(!(param %in% names(dots) | paste0(param, "_abs") %in% names(dots))) {
            stop("Require argument: ",param," = ... or ",param,"_abs = ...")}
        if(param %in% names(dots)) {
            par_val <- x$n * dots[[param]]
            } else {
            par_val <- dots[[paste0(param,"_abs")]]
        }
        if(!(is.numeric(par_val) & (length(par_val) == 1 | length(par_val) == x$k))) {
            stop(param, " must be a numeric literal or k-length vector")
        }
        if(length(par_val)==1) par_val <- par_val * x$pi

        # Warn user in case accidentally specified (%) as (abs) or (abs) as (%).
        # If user has entered Infinity, we can assume this is not by accident.
        if(is.null(dots[["silent"]]) || dots[["silent"]] == FALSE) {
            if(max(par_val) >= 10^5 & is.finite(max(par_val))) {
                message(param," values above 10^5 - extremely strong prior")}
            if(max(par_val) <= 3 ) {
                message(param," values below 3 - extremely weak prior")}
        }
        
        # Ensuring hyperparameters obey distribution constraints on Dirichlet / Inverse Wishart.
        if(param == "alpha" & min(par_val) < 1) {
            message("alpha violates prior distribution constraints - values have been adjusted.")
            par_val <- par_val/min(par_val)}
        if(param %in% c("beta", "nu") & min(par_val) < x$d) {
            message("beta/nu violates prior distribution constraints - values have been adjusted.")
            par_val <- par_val + x$d - min(par_val)}

        return(par_val)
    }
    
    # PROCESS ALL PARAMETERS 
    alpha <- processParam(dots, "alpha")
    if(length(grep("beta",names(dots))) > 0) {
        beta <- processParam(dots, "beta")
    } else {
        beta <- alpha
    }
    if(length(grep("nu",names(dots))) > 0) {
        nu <- processParam(dots, "nu")
    } else {
        nu <- beta
    }
    
    # Denote infinite values in table
    infinities <- data.frame(alpha = numeric(x$k), beta = numeric(x$k), nu = numeric(x$k))
    infinities$alpha <- is.infinite(alpha)
    infinities$beta <- is.infinite(beta)
    infinities$nu <- is.infinite(nu)
    for(cI in 1:3) {
        if(any(infinities[ ,cI] == TRUE) & any(infinities[ ,cI] == FALSE)) {
            stop("Infinite hypervalues not consistent: ",colnames(infinities)[cI])
        }
    }

    return(list(alpha=alpha, beta=beta, nu=nu, infinities=infinities))
}


