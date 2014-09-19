plot.clsMR <-
function(x, ...) {
    Call <- match.call(expand.dots = TRUE)
    stopifnot(class(x) == "clsMR")
    
    #If user specifies clusters from kept cluster data:
    if(length(Call) >= 3) {
        y <- as.cTurn(x, ...)
        plot(y)
    } else {
    
    #Else print agglomeration schedule in igraph
    	require(igraph)
        g <- graph.edgelist(el = na.omit(cbind(1:(nrow(x$agglom)-1), x$agglom[-nrow(x$agglom), "parent"])))
    	if(nrow(x$agglom) - length(V(g)) > 0) g <- add.vertices(g, nrow(x$agglom) - length(V(g)))
        v.size <- (x$agglom[ ,"n"]/max(x$agglom[ ,"n"]))^0.8
        opar <- par()$mar
        par(mar = c(0,0,0,0))
        plot(g, vertex.size = 28*v.size,layout = x$agglom[ ,c("r", "position")], vertex.color = "darkseagreen",
             vertex.frame.color = "grey45", edge.color = "palegreen3")
        par(mar = opar)
    }
}
