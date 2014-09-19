summary.cTurn <-
function(object, ...) {
    stopifnot(class(object) == "cTurn")
    print(object)
    if(object$n == 0) return()
    cat("Cluster size distribution summary: \n\n", sep="")
    print(summary(object$table$Freq))
    cat("\n")
    tbl <- t(matrix(as.numeric(as.matrix(object$table)), ncol = 2))
    n.total <- length(object$cluster)
    
    tbl <- rbind(tbl, round(100*tbl[2, ]/n.total,1))
    colnames(tbl) <- object$table$res    
    tbl <- tbl[2:3,,drop = F]
    rownames(tbl) <- c("n   ","%   ")
    print(tbl)
    
    if("cluster" %in% names(object)) {
        cls.means <- object$means
        
        colnames(cls.means) <- paste0("dim", 1:ncol(cls.means))
        rownames(cls.means) <- 1:object$k
        cat("\nCluster Means:\n")
        print(t(cls.means))
    }
}
