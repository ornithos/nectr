.nectr.getData <- function(x, append.cls = FALSE, omit.na = FALSE) {
    if(class(x) == "cTurn") {
        dat <- .nectr.fetchDataset(x$dataset)
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
            message(".nectr.getData: append.cls and omit.na are unavailable when x is not a cTurn object.")}
        return(.nectr.fetchDataset(x))
    } else if(class(x) == "clsMR") {
        if(append.cls | omit.na) {
            message(".nectr.getData: append.cls and omit.na are unavailable when x is not a cTurn object.")}
        return(.nectr.fetchDataset(x$dataset))
    } else {
        stop("no method for class '", class(x), "'.")
    }
}

.nectr.fetchDataset <- function(x) {
    dat <- try(get(x, envir = .GlobalEnv), silent=TRUE)
    if(inherits(dat, "try-error")) stop("Dataset '", x,"' does not exist in .GlobalEnv")
    return(dat)
}

.nectr.checkParams <- function(x) {
    slots <- names(x)
    if(class(x) == ".nectr.suffStats") { 
        err.msg <- "Object is not of valid GMM parameter type"
        if(!(length(slots) == 7 & length(union(slots, c("pi","mu","S","k","d","n","dsn")))==7)) {
            stop(err.msg, ": (incorrect object structure)")}
        if(any(x$pi < 0)) stop(err.msg, ": (all elements of pi should be nonnegative)")
        if(sum(x$pi) > 1 + 10^-9 | sum(x$pi) < 1 - 10^-9) stop(err.msg, ": (pi doesn't sum to 1)")
        if(!(length(x$k) == 1 & length(x$d) == 1 & is.numeric(x$k) & is.numeric(x$d))) {
            stop(err.msg, ": (missing num clusters k or dimension d)")}
        if(!(class(x$mu) %in% c("matrix","data.frame"))) {
            stop(err.msg, ": (mu must be specified as columnwise ", x$k, " by ", x$d, " matrix.)")}
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
    } else if(class(x) == ".nectr.Hyper") {
        err.msg <- "Object is not a valid hyperparameter object"
        slots <- names(x)
        if(!(length(slots) == 4 & length(union(slots, c("alpha","beta","nu","infinities")))==4)) {
            stop(err.msg, ": (incorrect object structure)")}
        l <- length(x$alpha)
        if(!is.data.frame(x$infinities)) stop(err.msg, ": (infinities must be a data.frame)")
        if(!(length(x$beta)==l & length(x$nu)==l & nrow(x$infinities)==l)) {
            stop(err.msg, ": (inconsistent dimensions)")}
        if(!(is.numeric(x$alpha) & is.numeric(x$beta) & is.numeric(x$nu))) {
            stop(err.msg, ": (alpha/beta/nu must be numeric vectors)")}
        if(ncol(x$infinities) != 3) stop(err.msg, ": (infinities should have 3 columns (alpha/beta/nu))")
        return(invisible(TRUE))
    } else {
        stop("Object is not of valid GMM parameter type: (incorrect class)")
    }
}

.nectr.checkDataset <- function(data) {
    # Check dataset is numeric matrix / dataframe
    dtype <- TRUE
    dftype <- TRUE
    if(class(data)=="matrix") {
        if(!(mode(data)=="numeric"|mode(data)=="logical")) {
            stop("Input data not numeric! \n\n")
        }
        dftype <- FALSE
    } else if(class(data)=="data.frame") {
        for(i in 1:ncol(data)) {
            if(!(mode(data[,i]) == "numeric" | mode(data[,i]) == "logical")) {
                stop(paste0("Column ",i," of data is not numeric! \n\n"))
            }
        }
    } else {
        dtype <- FALSE
        dftype <- FALSE
    }
    # TRUE = data.frame
    return(c(dtype, dftype))

}

# User-friendly function to generate the hyperparameters. Can be expressed as only 2 literals,
# or with alpha, beta, nu specified individually for each k. Handles both (%) and (abs#)
.nectr.getHyper <- function(x, ...) {
    dots <- list(...)
    if(class(x) != ".nectr.suffStats") stop("x must be of type .nectr.suffStats")
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
            if(max(par_val) >= 5*10^5 & is.finite(max(par_val))) {
                message(param," values above 5e+5 - strong prior")}
            if(max(par_val) <= 3 ) {
                message(param," values below 3 - very weak prior")}
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

    out <- list(alpha=alpha, beta=beta, nu=nu, infinities=infinities)
    class(out) <- ".nectr.Hyper"
    return(out)
}



as.SuffStats <- function(x) {
    slots <- names(x)
    if(!(length(intersect(slots, c("pi","mu","S","k","d","n","dsn")))==7)) {
        stop("Object does not contain all the required elements for SuffStats class.")}
    for(nm in slots) {
        if(nm %in% c("pi","mu","S","k","d","n","dsn")) next
        x[nm] <- NULL
    }
    class(x) <- ".nectr.suffStats"
    return(x)
}

