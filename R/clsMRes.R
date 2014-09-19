clsMRes <-
function(data, keep, r.start = NA, r.max = Inf, ...) {
    
    mc <- match.call()
    dataset.name <- toString(as.list(mc)$data)
    
    #***********************************
    #  INITIALISATION
    #***********************************
    
    r <- 1        #Initialise resolution
    i <- 1        #Initialise iteration number
    orig.n <- nrow(data)
    p <- ncol(data)
    
    #Generic results printer in "pretty" format
    print.summary <- function(x) {
        print.out <- c(formatC(x[1], digits = 3, width = 6, format = "f"),
                       formatC(x[2], digits = 0, width = 6, format = "f"),
                       paste0(formatC(100*x[3]/orig.n, digits = 2, width = 6, format = "f"),"%"))
        names(print.out) <- NULL
        print(noquote(print.out))
    }
    
    #Table "headers" for summary output
    print(noquote(formatC(c("r","k","n"), width = 6, format = "s")))
    
    if(is.na(r.start)) {
        #***********************************
        #  FIND S(Inf)
        #***********************************
        
        #~~~~~~RESOLUTION DESCENT SEQUENCE~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        base.seq <- c(seq(0.8, 0.2, by = -0.2), 1/seq(10, 50, by = 10), 1/seq(80, 250, by = 50))
        dataset.sd <- sd(as.numeric(as.matrix(data)))
        crude.guess <- c(0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000)[findInterval(dataset.sd, 10^seq(-2,3))]
        #-----TEST GUESS----------------------------------------------
        #--If our initial guess is already at S(inf) we don't want to climb up linearly from
        #--S(Inf) if our resolution starts way too low.
        while(TRUE) {
            curr <- clsTurnRes(data, crude.guess, summarise = 2, ...)
            if(!(curr$stats[3] <= 1 & curr$stats[2] <= 0.01*orig.n)) break
            crude.guess <- crude.guess*1.5
        }    
        descent.seq <- crude.guess*base.seq
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        prev <- c(NA, NA)
        
        #Descend until S(Inf).
        for(r in descent.seq) {
            curr <- clsTurnRes(data, r, summarise = 2, ...)
            print.summary(c(curr$stats["r"], curr$stats["k"], curr$stats["n"]))
            if (curr$stats["n"] <= 0.04*orig.n | 
                    (curr$stats["k"] <= 4 & curr$stats["n"] <= 0.1*orig.n)) break
            if(!is.na(prev[2])) {
                if((curr$stats["n"] - prev[2])/orig.n > -0.01) break
            }
            prev[2] <- prev[1]
            prev[1] <- curr$stats["n"]
        }
        
        cat(noquote(sprintf("Found S(Inf) at r = %#.4f", r)), "\n\n")
    } else {
        r <- r.start
    }
    
    #***********************************
    #  Initialise for output
    #***********************************
    
    #~~~~~~Set up output data.frames~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nodes.max <- 2000
    nodes.curr <- 1
    nodes <- matrix(NA, nrow = nodes.max, ncol = 4)
    colnames(nodes) <- c("r","n","parent","branches")
    
    turn.curr <- 1
    turn <- matrix(NA, nrow = 250, ncol = 7)
    colnames(turn) <- c("r", "n", "k", "t", "tm", "na", "tna")
    if(keep) cls.keep <- matrix(NA, nrow = orig.n, 10)
    all.means <- matrix(NA, 0, p)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    #***********************************
    #  FIND S(1)
    #***********************************
    
    #~~~~~~RESOLUTION ASCENT SEQUENCE~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ascent.seq <- r*(1 + seq(0.7,70, by = 0.7))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    #Ascend until S(1).
    for(r in ascent.seq) {
        curr <- clsTurnRes(data, r, summarise = 2, ...)
        print.summary(c(curr$stats["r"], curr$stats["k"], curr$stats["n"]))
        k <- curr$stats["k"]
        
        
        #Add to TURN results matrix
        turn[turn.curr, ] <- curr$stats
        turn.curr = turn.curr + 1
        
        #NODES/EDGES (if k == 0, there is nothing to add)
        if (k > 0) {
            #Expand output datasets if necessary.
            if(nodes.curr + k > nodes.max) {
                nodes <- rbind(nodes, matrix(NA, nrow = 1000, ncol = 5))
            }
            
            new.node.ids <- (nodes.curr):(nodes.curr + k-1)
            #Add node data into table
            nodes[new.node.ids, ] <- cbind(r, curr$freq, rep(NA,k), rep(0,k))
            
            #Add Edges from Previous Iteration in agglom schedule
            if (nodes.curr > 1) {
                prev.nodes <- (prev.curr):(nodes.curr - 1)
                movement <- table(prev.cls, curr$cluster)
                cls.into <- apply(movement, 1, which.max)
                
                e <- try(nodes[prev.nodes, "parent"] <- cls.into + nodes.curr-1, silent = TRUE)
                if(inherits(e, "try-error")) browser()
                
                branch.lkp <- as.data.frame(table(nodes[prev.nodes, "parent"]))
                branch.new <- branch.lkp[match(new.node.ids, branch.lkp[ ,1]), 2]
                nodes[new.node.ids, "branches"] <- ifelse(is.na(branch.new),0,branch.new)
                prev.curr <- nodes.curr
                nodes.curr <- nodes.curr + k
                
            } else {
                
                #First iteration only
                prev.curr <- nodes.curr
                nodes.curr <- nodes.curr + k
            }
            
            if(keep) {
                if((turn.curr-1) > ncol(cls.keep)) {
                    cls.keep <- cbind(cls.keep, matrix(NA, orig.n, 5))
                }
                cls.keep[ ,turn.curr-1] <- curr$cluster
            }
            all.means <- rbind(all.means, curr$means)
            prev.cls <- curr$cluster
        }
        
        #Exit criterion.
        if ((k == 1 & curr$stats["n"] >= 0.85*orig.n)|
                (curr$stats["n"] > 0.92*orig.n)|
                (r >= r.max)) {
            nodes <- nodes[1:(nodes.curr-1),]
            turn <- turn[1:(turn.curr-1),]
            break
        }
    }
    
    msg <- sprintf("Found S(1) at r = %#.4f. %0.2f%% of dataset clustered", r, 100*curr$stats["n"]/orig.n)
    cat(noquote(msg))
    
    if(keep) cls.keep <- cls.keep[ ,1:(turn.curr-1)]
    
    #***********************************
    #  GENERATE AGGLOMERATION SCHEDULE
    #***********************************
    
    rs <- rev(turn[ ,"r"])
    nodes <- cbind(nodes, position = NA)
    
    #Iterate from top of tree to leaves
    for(i in 1:length(rs)) {
        
        r <- rs[i]
        nodes.sel <- which(nodes[ ,"r"] == r)
        
        #If top of the tree (and only one node), then can shortcut to 0.5.
        if(i == 1 & length(nodes.sel) == 1) {
            nodes[nodes.sel, "position"] <- 0.5
            next
        }
        
        sum.branches <- sum(pmax(nodes[nodes.sel, "branches"], 1))
        
        #Get order of level above
        if(i > 1) {
            nx.nodes <- which(nodes[ ,"r"] == rs[i-1])
            nx.order <- order(nodes[nx.nodes, "position"])
            above.order <- match(nodes[nodes.sel, "parent"], nx.nodes[nx.order])
            above.order <- ifelse(is.na(above.order), length(nx.nodes) + 1, above.order)
        } else {
            above.order <- 1:length(nodes.sel)
        }
        
        #Order first by number of branches, and subsequently, by order at next level.
        #nodes.order <- nodes[order(tree$branches[nodes])]
        #nodes.order <- nodes[order(above.order)]
        
        #**** Re-order within each above.order - we want nodes with most branches to be higher
        #***  if less than mid point, and lower if greater than mid point
        nx.num <- length(unique(above.order))
        for(a in 1:nx.num) {
            i <- unique(above.order)[a]
            sub.select <- which(above.order == i)
            decr = T
            if(a <= nx.num/2) decr = F
            sub.order <- order(nodes[nodes.sel[sub.select], "branches"], decreasing = decr)
            nodes.sel[sub.select] <- nodes.sel[sub.select[sub.order]]
        }
        nodes.order <- nodes.sel[order(above.order)]
        
        #Calculate position
        p <- cumsum(pmax(nodes[nodes.order, "branches"], 1))
        nodes[nodes.order, "position"] <- p/(sum.branches + 1) 
    }
    
    out <- list()
    out$agglom <- nodes
    out$turndata <- turn
    out$dataset.name <- dataset.name
    out$means <- all.means
    out$res.step <- ascent.seq[2] - ascent.seq[1]
    out$n <- orig.n
    if(keep) out$keep <- cls.keep
    class(out) <- "clsMR"
    return(out)
}



"[.clsMR" <- function(x, y, ...) {
	x$agglom[y, c("r","n")]
	}