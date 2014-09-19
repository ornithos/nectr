\name{getEMGPs}
\alias{getEMGPs}
\alias{getEMClusters}

\title{Expectation - Maximisation Algorithm for Gaussian Mixtures}
\description{
	Implementation of Expectation - Maximisation (EM) Algorithm. Many implementations exist already within R. The
	reason that this package includes another is for scalability reasons and also some modifications. It is intended
	to be used in series with TURN-RES, and the analyst can lock the clusters found through this in place.
}

\usage{
getEMGPs(data, centers, cls.prob = NULL, move.restrict = NA, eps.target = 0.1,
					trunc = 0.5, max.iter = 100, silent = FALSE)
				
getEMClusters(data, mu = NA, sigma = NA, phi = NA, params = NA)
}
\arguments{
  \item{data}{\emph{required}. A numeric data frame or matrix where each column is a dimension to be clustered over.}
  \item{centers}{\emph{required}. List of vectors, where each vector is an initialisation for a cluster center. Implicit in
			this argument is the number of clusters.}
  \item{cls.prob}{Optional initialisation of cluster prior probabilities. For k clusters, defaults to n/k}
  \item{move.restrict}{Optional restriction of cluster center update. Clusters are limited to the nearest possible point 
			during the update process \emph{within a locus of} \code{move.restrict} around the center initialisation}
  \item{eps.target}{Optional convergence criteria. Algorithm has converged when likelihood function changes < \code{eps.target}}
  \item{trunc}{Optional paramater to remove outliers from covariance estimation. Datapoints with < trunc/n probabilities are
			removed from the calculation.}
  \item{max.iter}{Optional iteration criteria. Maximumm number iterations allowed to converge before algorithm exit.}
  \item{silent}{Optional parameter. By default, algorithm gives console update after every iteration.}
}

\details{No default initialisation of gaussian means is implemented, neither are the number of clusters estimated. The
		implementation is built to work alongside the exploratory clustering available through TURN-RES \link{clsTurnRes} and so
		it is assumed that the user has determined both the number of clusters and their centers already.
		
		The major modification to the algorithm is allowing the user to restrict how far the clusters can be moved from their
		initial center. A value of \code{move.restrict = 0} locks the centers in place, and a value of \code{move.restrict = 5}
		prevents the algorithm from moving the cluster centers further than a (euclidean) distance of 5 from their initial
		specification. Note that the convergence behaviour has not been analysed with this modification and the author makes no
		guarantees of convergence. In practice however, all known attempts have converged, but monotonicity of optimisation is
		violated.
		
		It has also been experimentally discovered that sufficient noise in the dataset can lead to pathological covariance matrices,
		which might have such a large volume that all noise is considered part of that cluster and the intended cluster is lost. The
		\code{trunc} argument has been implemented in order to 'trim' the covariance matrix input such that outliers are removed. The
		precise implementation is that any datapoint that has a probability of < trunc/n of belonging to any of the k clusters is
		removed from the covariance calculation of that iteration. Again, no guarantees are therefore made of convergence.
		}

\value{
	A list with slots for the cluster parameters (\code{phi},\code{mu},\code{sigma}), the \code{cluster} assignment vector, and the
	\code{cluster.prob} vector. The cluster.prob vector is the density of the assigned cluster at the given datapoint; this can be used
	as a measure of strength of cluster membership. The list is unclassed; it is a generic object to be used as desired by the end user;
	no further functionality is given.
}

\note{
	This implementation is heavily indebted to the package \code{\link[mvtnorm:dmvnorm]{mvtnorm}} for providing a fast calculation of multivariate 
	gaussian density.
}


\author{
	Alex Bird, \email{alex.bird@boots.co.uk}
}

\seealso{\code{\link{getEMClusters}} for scoring the data into clusters using the parameters estimated from this function.
		\code{\link{formaliseClusters}} for guided use of these functions or if using input from \code{\link{clsMRes}} /
		\code{\link{clsTurnRes}} objects.}

\examples{
x = rnorm(30,0,1);
mean(x) #approx 0.5. Nothing to do with this!
}
