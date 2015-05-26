as.cTurn <- function(x, ...) {
    if(class(x) == "cTurn") return()
    if(class(x) == "clsMR") {
        if(is.null(x$keep)) stop("Mres object must have keep = TRUE to perform operation")
        Call <- match.call(expand.dots = TRUE)
        if(length(Call) < 3) stop("Specify nodes to collapse into cTurn (eg. as.cTurn(x, c(1,11,21))")
        nodes <- eval(Call[[3]])
        if(any(!(unlist(nodes) %in% 1:nrow(x$agglom)))) stop("nodes must be in range 1:", nrow(x$agglom))
        nodes <- as.list(nodes)
        
        y <- list()
        y$k <- length(nodes)
        y$dataset <- x$dataset
        
        #Get column/resolution of each cluster specified
        cumcls <- cumsum(x$turndata[ ,"k"])
        run.num <- list()
        for(a in 1:length(nodes)) {
            for(b in nodes[[a]]) {
                for(s in 1:length(cumcls)) {
                    if(b <= cumcls[s]) {
                        if(length(run.num) < a) {
                            run.num[a] <- s
                        } else {
                            run.num[[a]] <- c(run.num[[a]], s)
                        }
                        break
                    }
                }
            }
        }
        
        ul.run.num <- unlist(run.num)
        ul.nodes <- data.frame(nodes = unlist(nodes), 
                               lkp = rep(1:length(nodes), sapply(nodes, function(z) length(z))))
        schedule.res <- sort(unique(ul.run.num))
        cumcls <- c(0, cumcls)
        cluster <- rep(NA,x$n)
        for(a in rev(schedule.res)) {
            nodes.x <- ul.nodes$nodes[which(ul.run.num == a)] - cumcls[a]
            cls.x <- x$keep[, a] %in% nodes.x
            cls.nodes <- x$keep[cls.x ,a] + cumcls[a]
            cls.nodes <- merge(data.frame(id = 1:length(cls.nodes), nodes = cls.nodes), ul.nodes)
            cluster[cls.x] <- cls.nodes[order(cls.nodes$id), 3]
        }
        y$cluster <- cluster
        y$table <- as.data.frame(table(cluster))
        names(y$table) <- c("res", "Freq")
        y$n <- sum(!is.na(y$cluster))
        y$r <- "mres"
        y$means <- byMeans(get(y$dataset, envir = .GlobalEnv), cluster, 1:length(nodes))
        class(y) <- "cTurn"
        return(y)
    }
    else if(class(x) == "kmeans") {
        Call <- match.call(expand.dots = TRUE)
        if(length(Call) < 3) stop("Requires string of dataset name (eg. as.cTurn(x, 'dataset')")
        dataset.name <- eval(Call[[3]])
        
        y <- list()
        y$k <- length(x$size)
        y$dataset <- dataset.name
        y$n <- length(x$cluster)
        y$cluster <- x$cluster
        y$table <- data.frame(res = 1:y$k, Freq = x$size)
        y$r <- "kmeans"
        y$means <- byMeans(get(dataset.name, envir = .GlobalEnv), cluster, 1:y$k)
        class(y) <- "cTurn"
        return(y)
    }
    else {
        stop("Unable to process object of class ", class(x), ". Requires clsMR or kmeans object")
    }
}


clsMerge <- function(x, node.list) {
    as.cTurn(x, node.list)
}
    