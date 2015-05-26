plot.cTurn <-
function(x, ...) {
    dots <- list(...)
    stopifnot(class(x) == "cTurn")
    
    #No additional arguments
    if(length(dots) == 0) {
        if(x$k > 49) stop("Will not attempt to plot more than 49 clusters")
        
        #NOTE: Fn uses original dataset - if you manipulate this after clustering, plot will be trashed.
        require(ggplot2)
    	require(reshape2)
    	nm <- x$dataset
        title <- paste("Parallel Co-ordinate Plot for Clusters in dataset:",nm)
        message("Plotting Parallel Co-ordinate plot - this may take some time...")
        
        #Scale down large datasets
        if(x$n > 2000) {
            scale <- ceiling(x$n/2000)
            cls.names <- as.numeric(as.character(x$table$res))
            for(i in 1:nrow(x$table)) {
                j <- cls.names[i]
                n.tf <- pmin(pmax(100,trunc(x$table$Freq[i]/scale)), x$table$Freq[i])
                rws <- which(x$cluster == j)
                rws <- rws[sample(x$table$Freq[i] ,n.tf)]
                if(i==1) {
                    plot.data <- cbind(.nectr.getData(nm)[rws, ], cls = rep(j, n.tf))
                } else {
                    plot.data <- rbind(plot.data, cbind(.nectr.getData(nm)[rws, ], cls = rep(j, n.tf)))
                }
            }
        } else {
            plot.data <- .nectr.getData(x)
            plot.data <- cbind(plot.data, cls = x$cluster[rws]) 
        }
        
        #Add IDs and put into "molten" state (see reshape2)
        if(is.matrix(plot.data)) plot.data <- as.data.frame(plot.data)
        data <- melt(cbind(ID=1:nrow(plot.data), plot.data), id.vars = c("ID","cls"), variable.name = "dimensions")
        
        #Perform plot
        plot.alpha <- 170/pmax(170,pmin(x$n,2000))
        p <- ggplot(data, aes(x = dimensions, y = value, group = ID, colour = factor(cls))) + geom_path(alpha = plot.alpha) 
        p <- p + facet_wrap( ~ cls) + theme(axis.text.x=element_text(angle = -90, hjust = 0.5, vjust = 0.5))
        print(p + ggtitle(title) + theme(legend.position="none"))

    }
    #Extra arguments - plot specific clusters
    else {
    
        if(length(dots) == 1) type <- "box" else type <- dots[[2]]
        i <- dots[[1]]
        if(!is.vector(i)) stop("Second argument in plot must be a vector")
        if(!(class(i) %in% c("numeric", "integer"))) stop("Second argument must specify cluster numbers")
        
        if (any(!(i %in% x$table$res))) stop(cat("Cluster ",i," does not exist in this object.\n", sep = ""))
        data <- .nectr.getData(x)
        
        require(ggplot2)
        require(reshape2)
        if(type == "box") {
            bpset <- melt(data[!is.na(x$cluster) & x$cluster == i, ])
            qplot(x = factor(bpset$variable), y = bpset$value, geom = "boxplot")
        } else if (type == "pcp") {
            rws <- !is.na(x$cluster) & x$cluster %in% i
            pcpset <- cbind(ID = 1:(sum(rws)), data[rws, ])
            if(nrow(pcpset)> 3000) pcpset <- pcpset[sample(nrow(pcpset), 3000), ]
            pcpset <- melt(pcpset, id.vars = "ID", variable.name = "dimensions")
            v.cls <- x$cluster[x$cluster %in% i]
            pcpset$cls <- v.cls[pcpset$ID]
            p <- ggplot(pcpset, aes(x = dimensions, y = value, group = ID, col = factor(cls))) + geom_path(alpha = 0.03) 
            p <- p + theme(axis.text.x=element_text(angle = -90, hjust = 0.5, vjust = 0.5))
            print(p + ggtitle(paste("Parallel Co-ordinates Plot of cluster(s)", paste0(i, collapse = ","))))
        } else {
            stop("Invalid chart type specified: ",type)
        }
    }
}