\name{nectr}
\alias{nectr}


\title{nectr: Non-Parametric Exploratory Clustering with Turn-Res}
\description{
	This package aims to facilitate exploratory clustering with large low-dimensional datasets. In this context, 'large' means 
	datasets of O(10^6) rows and 5-10 columns. There are many clustering packages available in R, but very few that cope with
	enterprise size data. This package is designed for the industry analyst; analysis is guided, interactive and as easy to use
	as possible. Two algorithms are implemented here: TURN-RES (Foss 2002), and EM clustering for Gaussian Mixtures.
}

\details{The following functions are available:
		\itemize{
		\item{\link{clsTurnRes}:	An implementation of TURN-RES. This is the workhorse algorithm of the package.}
		\item{\link{clsMRes}:		A function which calls \code{clsTurnRes} over increasing resolution in order to determine
									the structure of the dataset. This alleviates the requirement to know parameters
									in advance.}
		\item{\link{formaliseClusters}:	A helper function which extracts clusters from the TURN-RES results and
									outputs the best fit Gaussian clusters. This is a guided use of the EM functions below.}
		\item{\link{getEMGPs}:			An implementation of EM clustering for Gaussian mixture models. This function
									returns merely the model parameters.}
		\item{\link{getEMClusters}:		Scores input data according to the GMM parameters obtained in \code{getEMGPs}.}
		}
		
		The TURN-RES algorithm is extremely fast (8-10 seconds for a 600k x 3 dataset on a 2.5 GHz i5 processor) and scalable: O(n log n). Density based in nature,
		the algorithm is not limited to convex or ellipsoidal clusters and noise below a threshold is excluded. A hierarchical view of the cluster 
		agglomeration in the algorithm is available by iterating through multiple values of the resolution parameter, which can be helpful in
		understanding the structure of the data.
		
		The TURN-RES algorithm suffers from a significant flaw. One often finds that clusters of interest are revealed across different values of the resolution
		parameter. The previously discovered clusters can be lost at a lower resolution, and amalgamating between resolutions is hard. A number of different
		approaches are available here, but within this package, a Gaussian Mixture Model has been used, which takes into account the changing density 
		throughout the cluster. This does of course reduce the flexibility in cluster shape somewhat - non-convex shapes in particular will be lost; however
		practically, this algorithmic choice has been helpful in other ways. Information regarding outliers is retained via the probability calculations and 
		unseen data can quickly and easily be scored - allowing clustering to be run on a subset of the original dataset if required.
		}


\examples{
\dontrun{
#Synthetic data
library(MASS)
mu <- list()
mu[[1]] <- c(4,6,3)
mu[[2]] <- c(1,0,-1)
mu[[3]] <- c(-4,-1,-2)
mu[[4]] <- c(2,-5,-4)
mu[[5]] <- c(-3,3,5)

sigma <- list()
sigma[[1]] <- matrix(c(2,0.2,0.8,0.2,1,0.1,0.8,0.1,1),3,3)
sigma[[2]] <- matrix(c(1,0,-1,0,2,0.4,-1,0.4,3),3,3)
sigma[[3]] <- matrix(c(3,0,0,0,0.3,0,0,0,2),3,3)
sigma[[4]] <- matrix(c(2,1.5,1,1.5,2,1,1,1,2),3,3)
sigma[[5]] <- diag(3)*1.8

data <- as.data.frame(rbind(mvrnorm(150000, mu[[1]], sigma[[1]]), mvrnorm(60000, mu[[2]], sigma[[2]]),
                            mvrnorm(52500, mu[[3]], sigma[[3]]), mvrnorm(135000, mu[[4]], sigma[[4]]), 
                            mvrnorm(105000,mu[[5]], sigma[[5]]),cbind(matrix(runif(75000*3,-10,10), 75000, 3))))
names(data) <- c("X", "Y", "Z")

#Run through Multi-Resolution TURN-RES
mres <- clsMRes(data, keep = TRUE)

#Plot agglomeration and Parallel Coordinate Charts
plot(mres)
plot(mres, c(4, 8, 9, 11))

#Obtain GMM parameters for centers desired
gmm.model <- formaliseClusters(mres, c(4, 8, 9, 11))
}
}
