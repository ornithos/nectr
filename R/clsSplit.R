clsSplit <-
function(x, cls.spec, k) {
    stopifnot(class(x) == "cTurn")
    if(!(is.vector(cls.spec) & class(cls.spec) == "numeric")) stop("cls.spec must be a numeric vector")
    if(length(k) > 1 || !is.vector(k)) stop("k must be a single value!")
    if(is.null(x$cluster)) stop("No cluster vector attached to cTurn object!")
    
    tmp.data <- get(x$dataset, envir = parent.env(environment()))[x$cluster %in% cls.spec, ]
    tmp.km <- kmeans(tmp.data, k)
    lkp <- c(cls.spec, seq(x$k+1, x$k+1+k-length(cls.spec), length.out = pmax(0,k-length(cls.spec))))
    x$k <- x$k - length(cls.spec) + k
    x$cluster[x$cluster %in% cls.spec] <- lkp[tmp.km$cluster]
    x$table <- as.data.frame(table(x$cluster))
    colnames(x$table)[1] <- "res"
    return(x)
}
