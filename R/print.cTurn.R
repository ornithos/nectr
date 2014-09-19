print.cTurn <-
function(x, ...) {
    stopifnot(class(x) == "cTurn")
    n.total <- length(x$cluster)
    n.pct <- 100*x$n/n.total
    cat("TURN-RES Clustering of '", x$dataset, "' dataset at r = ", x$r, "...\n\n", sep="")
    cat(paste(sprintf("Clustered %.1f%% of the dataset; n = ", n.pct),format(n.total, big.mark = ","),
              ", k =",x$k, "\n"))
}
