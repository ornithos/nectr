

# NOTE THAT ALPHA, ALPHA_ABS, BETA, BETA_ABS ARE IMPLICITLY IN HERE TOO
clsGMM <- function(x, clusters=NULL, noise=NULL, add.pi = NULL, add.mu = NULL, add.cov = NULL, ...) {
    
    args <- match.call(expand.dots = TRUE)
    if(!(class(x) %in% c("cTurn", "clsMR"))) 
        stop("object class '", class(x),"' not recognised. Supply either cTurn or clsMR object.")
        
    SS <- .nectr.getSuffStats(x, clusters)
    hyp <- .nectr.getHyper(SS, ...)
    
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



# Functionality:
# - add noise cluster.. ie large variance component, with mean as the mean of the data.
#
# - how do we play good access to hyperparameters? I think adding in components as well
#   makes everything too confusing. Suggest putting the augmented components - including
#   noise as a separate stage. "SPECIFY MODEL" stage. clsSpecifyModel?
#
# - And remember ease of use for hyperparams - ideally specify one or at most two numbers
#   perhaps 0 = no optimisation and 1 = full optimisation (ie use priors only as initialisation
#   save us from going for 1 -> Infty.
#
# - Think about external access to GMM function. To be fair this could be done via
#   Specify Model stage.
#
# - Final other access to GMM algorithm will be via predict.nectrGMM.