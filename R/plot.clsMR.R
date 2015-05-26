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
        g <- graph.edgelist(el = na.omit(cbind(1:(nrow(x$agglom)-1), x$agglom[-nrow(x$agglom), "child"])))
    	if(nrow(x$agglom) - length(V(g)) > 0) g <- add.vertices(g, nrow(x$agglom) - length(V(g)))
        v.size <- (x$agglom[ ,"n"]/max(x$agglom[ ,"n"]))^0.8
        opar <- par()$mar
        par(mar = c(0,0,1.5,0))
        plot(g, vertex.size = 28*v.size,layout = x$agglom[ ,c("r", "position")], vertex.color = "#BE4C24",
             vertex.frame.color = "#CDB35D", edge.color = "#A68E86", vertex.label.cex=1.4, main =
                 "Cluster agglomeration schedule")
        par(mar = opar)
    }
}
