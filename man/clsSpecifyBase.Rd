\name{clsSpecifyBase}
\alias{clsSpecifyBase}

\title{Specify New TURN-RES Base by Agglomerating Existing Clusters}
\description{
	A helper function which creates a base for a new TURN-RES iteration. This function will agglomerate
	specified clusters and return a cluster membership vector that can be used in the \code{base.cls} argument
	in \link{clsTurnRes}.
}

\usage{
clsSpecifyBase(x, clusters, trim = 0.05)
}
\arguments{
  \item{x}{a \link{cTurn} object containing the clustering to be used.}
  \item{clusters}{a vector the same length as the number of clusters \code{k} of \code{x} which acts as an agglomeration
		map. For example, if x has 5 clusters, the vector \code{c(1,2,1,3,3)} would agglomerate both the first and third 
		cluster and the fourth and fifth cluster, leaving the second as is.}
  \item{trim}{\emph{optional} Specifies whether the existing clusters are to be trimmed. By default, this value is set to
		0.05, which will trim 0.025 quantiles off either end of each cluster in every dimension. Clearly setting this 
		parameter to 0 will remove this feature.
	}
}
\details{This is designed to be used with the \code{base.cls} argument of a \link{cTurn} object. Hence trimming each cluster
		can be helpful if some of the edge datapoints are not representative of the rest of the cluster.}

\value{
	A vector corresponding to a grouped and trimmed version of the cluster membership vector of a \code{cTurn} object.
}

\author{
	Alex Bird, \email{alex.bird@boots.co.uk}
}

\examples{
x = rnorm(30,0,1);
mean(x) #approx 0.5. Nothing to do with this!
}
