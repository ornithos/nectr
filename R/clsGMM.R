clsSpecifyModel <- function(x = NULL, clusters=NULL, noise.pct=0, add.pi = NULL, 
                            add.mu = NULL, add.cov = NULL, dataset = NULL) {
    
    cTurn <- FALSE
    add <- FALSE
    if(length(noise.pct) != 1) stop("noise.pct must be a single numeric value between 0 and 1")
    if(noise.pct<0 || noise.pct>1) stop("noise.pct must be a value between 0 and 1")
    if(!is.null(x)) {
        cTurn <- TRUE
        if(!(class(x) %in% c("cTurn", "clsMR", ".nectr.suffStats"))) 
                stop("object class '", class(x),"' not recognised. Supply either cTurn or clsMR object.")
        if(class(x) != ".nectr.suffStats")
            ct <- x
            x <- .nectr.getSuffStats(x, clusters)
        # x becomes sufficient statistics - save the input object as ct
    }
    if(!(is.null(add.pi) && is.null(add.mu) && is.null(add.cov))) add <- TRUE
    if(!cTurn && is.null(dataset)) 
        stop("If no cTurn/clsMR object specified, dataset reference must be given.")
    if(!cTurn && !add)
        stop("What do you think you're trying to accomplish here?")
    
    # Adding in extra components to model if applicable.
    if(add) {
        if(cTurn) {
            x <- list()
            if(noise.pct!=0) {
                x$pi <- c(add.pi, noise.pct)
            } else {
                x$pi <- add.pi
            }
            if(sum(x$pi) > 1) {
                x$pi <- x$pi/sum(x$pi)
            } else if(sum(x$pi) < 1) {
                warning("Probabilities (pi) and noise.pct sum to less than 1. Rescaled.")
                x$pi <- x$pi/sum(x$pi)
            }
            x$mu <- add.mu
            x$S <- add.cov
            x$k <- length(add.pi)
            x$d <- ncol(dataset)
            x$n <- nrow(dataset)
            x$dsn <- toString(as.list(match.call())$dataset)
            check <- .nectr.checkDataset(dataset)
            if(!check[1]) stop("Dataset specified is not a numeric matrix or data.frame")
            if(!check[2]) dataset <- as.data.frame(dataset)
            class(x) <- ".nectr.suffStats"
            .nectr.checkParams(x)
        } else {
            add.pi.sum <- noise.pct + sum(add.pi)
            if(add.pi.sum  >= 1) 
                stop("Noise and add.pi probabilities must sum to < 1")
            x$pi <- c(x$pi*(1-add.pi.sum), add.pi)
            if(noise.pct>0) x$pi <- c(x$pi, noise.pct)
            
            k <- length(add.pi)
            if(k ==1) {
                if(!((is.vector(add.mu) && length(add.mu) == x$d) || 
                         (is.matrix(add.mu) && prod(dim(add.mu))==x$d)))
                    stop("mu must be a vector of length ",x$d)
            } else {
                if(!((is.matrix(add.mu)||is.data.frame(add.mu)) && ncol(add.mu)== x$d && 
                         nrow(add.mu)==k))
                    stop("mu must be a matrix with ",k, " rows and ",x$d, " columns.")
            }
            x$mu <- rbind(x$mu, as.vector(add.mu))
            x$S <- c(x$S, add.cov)
            x$k <- x$k + k
            .nectr.checkParams(x)
        }
    }
        
    # Adding in noise cluster
    if(noise.pct > 0) {
        x$mu <- rbind(x$mu, apply(x$mu, 2, median))
        if(cTurn)
            dataset <- .nectr.getData(ct)
        x$S <- c(x$S, list(cov(dataset)))
        if(!add) {
            # Note that if components were added, the noise probabilities were
            # updated already (for convenience).
            x$pi <- c(x$pi*(1-noise.pct), noise.pct)
        }
        x$k <- x$k + 1
    }
    .nectr.checkParams(x)
    
    # Finish up object
    if(x$k == 1 && noise.pct == 0) stop("Model is redundant: only one cluster specified.")
    x$cluster <- NA
    x$noise <- ifelse(noise.pct>0, x$k, NA)
    x$llh_ascent <- NA
    x$monotone <- NA
    x$iter <- NA
    x$fit <- FALSE
    class(x) <- ".nectr.GMM"  
    return(x)
}


# NOTE THAT ALPHA, ALPHA_ABS, BETA, BETA_ABS ARE IMPLICITLY IN HERE TOO
clsGMM <- function(x,  ...) {
    
    args <- match.call(expand.dots = TRUE)
    if(!(class(x) == ".nectr.GMM")) 
        stop("object x must be an unfitted .nectr.GMM object.\n",
             "These can be created using the clsSpecifyModel function.")
        
    SS <- as.SuffStats(x)
    hyp <- .nectr.getHyper(SS, ...)
    
    # Remove unnecessary arguments from the input to fitting function
    dots <- as.list(substitute(list(...)))[-1]
    reqd <- as.list(args(.nectr.fitGMM))[-1]
    keep <- names(dots) %in% names(reqd)
    dots <- dots[keep]
    dots$params <- SS
    dots$hyperparams <- hyp
    
    if(!is.null(silent)) cat("Performing EM iterations to fit specified GMM:\n")
    fit <- do.call(".nectr.fitGMM", dots)
    return(fit)
}



# Functionality:
# - need to add noise cluster argument into .nectr.fitGMM.
#   The hope being that the cluster output also gives the next nearest 'non-noise' cluster
#   should the user wish to classify by this. Maybe an addl column onto the cluster df
#   perhaps also add probability again onto here... since gen model (in theory) we can use bayes
#   to give 'probability of noise'..!

# - Final other access to GMM algorithm will be via predict.nectrGMM.

# - clsSplit is an interesting functionality to add via kmeans

# - THEN... need to add summary and plot methods for .nectr.fitGMM