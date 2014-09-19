# dyn.load("./src/nb_all_in1.dll")
# dyn.load("./src/adjTreeCpp.dll")
# dyn.load("./src/clsMeans.dll")

clsTurnRes <-
function(data, r, summarise = F, min.size = "Auto", base.cls = "None", phi = 0.8) {
    
    #***************************************************
    #    DATA PREP
    #***************************************************
    
    #Catch input errors
    if(class(data)=="matrix") {
        if(!(mode(data)=="numeric"|mode(data)=="logical")) {
            stop("Input data not numeric! \n\n")
        }
        message("Converting input data to dataframe..")
        data <- as.data.frame(data)
    } else if(class(data)=="data.frame") {
        for(i in 1:ncol(data)) {
            if(!(mode(data[,i]) == "numeric" | mode(data[,i]) == "logical")) {
                stop(paste0("Column ",i," is not numeric! \n\n"))
            }
        }
    } else if(class(data) != "cTurn") {
        stop("data is not a matrix, dataframe or cTurn object! \n\n")
    }
    
    n <- ifelse(class(data) == "cTurn", length(data$cluster), nrow(data))
    
    if(class(base.cls) == "character") {
        if(base.cls != "None") stop("When specifying base clusters, must be in vector format!\n\n")
    } else if (!(class(base.cls) %in% c("integer","numeric"))) {
        stop("Base clusters must be specified in vector format!\n\n")
    } else if (length(base.cls) != n) {
        if((class(data) == "cTurn") &
               (length(base.cls) < pmin(pmax(10, n/100),100))) {
            if(!(all(base.cls == unique(base.cls)) & all(base.cls %in% data$table$res))) {
                stop("Base vector must either be same length as dataset, or be specifying valid cluster numbers")
            }
            #Convert base cluster specifcation to base cluster vector
            base.cls <- ifelse(data$cluster %in% base.cls, data$cluster, 0)
        } else {
            stop("Base clusters vector must be same length as data!\n\n")
        }
    }
    
    #If passed a cTurn object, need to extract dataset. Else record name from arguments
    if(class(data) == "cTurn") {
        dataset.name <- data$dataset
        data <- get(dataset.name, envir = parent.env(environment()))
    } else {
        #Get data object name
        mc <- match.call()
        dataset.name <- toString(as.list(mc)$data)
    }
    
    # Function Begin -------------------------------------------
    
    #Append IDs
    data <- cbind(ID = 1:n, data)
    
    #Do Resolution scaling
    data[ ,-1] <- trunc(data[ ,-1]/r)
    
    #Initialise vars
    p <- ncol(data) - 1 #(Ignore ID column)
    if(is.character(min.size)) {
        if(min.size == "Auto") noise.size <- max(round(n/100), 1)
    } else if(is.numeric(min.size)) {
        noise.size <- min.size - 1
    } else {
        stop("min.size must be numeric or 'Auto'.\n")
    }
    
    
    #Order dataset cyclically - each column from last to first
    permutation <- 1:n
    for(i in 1:p) {
        permutation <- permutation[order(data[permutation ,i+1])]
    }
    data <- data[permutation, ]
    
    #***************************************************
    #    CALCULATE NEIGHBOURS ####
    #***************************************************
    
    #Calculate neighbours - left and right in all dimensions and find density and closeness.
    #  computation done in C++ (~100x faster than R).
    nb.c <- .C("calculateNeighbours", ID = as.integer(data$ID), spatial = as.integer(as.matrix(data[ ,-1])), 
               nbs = as.integer(rep(0L, n*2*p)), SS = as.double(rep(0, n)), type = as.integer(rep(0,n)),
               n = as.integer(n), p = as.integer(p), phi = as.double(phi))  
    
    
    #***************************************************
    #   COMPUTE CLUSTERS                            ####
    #***************************************************
    
    #Housekeeping and variable initalisation for C procedure
    type <- nb.c$type
    nb <- matrix(nb.c$nbs, n, 2*p)
    if(is.character(base.cls)) {
        base.cls.specd <- FALSE
        base.cls <- rep(0, n)
        cls.order <- sample(n,n) - 1
    } else {
        base.cls.specd <- TRUE
        base.cls <- ifelse(is.na(base.cls), 0, base.cls)
        base.cls.size <- order(order(as.data.frame(table(base.cls))$Freq))
        base.cls.size <- c(length(base.cls.size)+1, base.cls.size)
        cls.order <- order(base.cls.size[base.cls+1])-1
    } 
    
    res <- .Call("doTurnTraversalCpp", nb, type, as.integer(cls.order), as.integer(base.cls))
    
    #SPECIFYING BASE CLUSTERS CAUSES HAVOC IN THE TREE TRAVERSAL ABOVE.-----------
    #If base clusters specified, merge all clusters described by base specification.
    #This bit of code could do with optimising when I have the time.
    if(base.cls.specd) {
        cls.sp.uq <- unique(base.cls)
        cls.sp.uq <- cls.sp.uq[cls.sp.uq!=0]
        mx.curr <- max(res)
        for(i in cls.sp.uq) {
            split.cls <- unique(res[base.cls == i])
            res[res %in% split.cls] <- mx.curr + i
        }
    }
    
    #Clearup.
    rm(cls.order, nb, base.cls)
    
    #***************************************************
    #    TIDY UP CLUSTERS AND COMPUTE STATS
    #***************************************************
    
    #Get list of cluster labels that have freq > noise.size.
    cls.freq <- as.data.frame(table(res))
    cls.freq$res <- as.numeric(levels(cls.freq$res)[as.numeric(cls.freq$res)]) #convert from factor
    noise <- cls.freq$res[cls.freq$Freq <= noise.size]
    cls.lkp <- cls.freq$res[!(cls.freq$res %in% noise)]
    
    #Lookup raw cls output - transform raw output to 1, .., k (NA if not in list above).
    if (length(cls.lkp) > 0) {
        res <- match(res, cls.lkp)
    } else {
        #No clusters large enough - all is noise.
        res <- rep(NA,n)
    }
    
    #Compute Stats
    results <- cbind(ID = nb.c$ID, tp = p/nb.c$SS, cluster = res)
    t.total <- sum(results[ ,2])
    
    results <- results[!(is.na(res)), , drop = F]
    if(length(dim(results)) == 2) {
        n.cls <- nrow(results)
    } else {
        n.cls <- 0
    }
    
    if(n.cls > 0) { 
        k <- length(cls.lkp) 
        cls.freq <- as.data.frame(table(res))
        #Calculate cluster means
        cls.means <- .Call("clsMeans", as.matrix(get(dataset.name, envir = parent.frame())), 
                           res, as.integer(cls.freq[,1]))
    } else { 
        k <- 0
        cls.freq <- data.frame(Var1 = 0, Freq = 0)
        cls.means <- matrix(NA,1,1)
    }    
    
    #******* OUTPUT *******************
    
    #TRUE == 1; 2 is the special case for function calls.
    if (summarise %in% 1:2) {
        
        
        gr <- list()
        gr$r <- r
        gr$n <- n.cls
        gr$k <- k
        gr$t <- sum(results[ ,2])
        if (gr$n == 0) { gr$tm <- 0 } else { gr$tm <- gr$t/gr$n }
        gr$na <- n - gr$n
        gr$tna <- t.total - gr$t
        freq <- cls.freq[ ,2]
        
        summ <- list()
        summ$stats <- unlist(gr)
        summ$freq <- freq
        summ$means <- cls.means
        
        if(summarise == 2) summ$cluster <- res
        return(summ)
        
    } else {
        
        full <- list()
        full$cluster <- res
        full$k <- k
        full$n <- n.cls
        full$table <- cls.freq
        full$means <- cls.means
        full$dataset <- dataset.name
        full$r <- r
        full$type <- type
        class(full) <- "cTurn"
        return(full)
    }
}

byMeans <- function(data, cluster.vector, unique.clusters = NA) {
    if(is.data.frame(data)) data <- as.matrix(data)
    if(!is.numeric(data)) stop("Data must be numeric")
    if(!is.vector(cluster.vector)) stop("cluster.vector must be a vector")
    if(class(cluster.vector) != "integer") stop("cluster.vector must be class integer")
    if(length(unique.clusters) && is.na(unique.clusters)) unique.clusters <- as.integer(na.omit(unique(cluster.vector)))
    return(.Call("clsMeans", data, as.integer(cluster.vector), as.integer(unique.clusters)))
}
