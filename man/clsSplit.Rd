\name{clsSplit}
\alias{clsSplit}

\title{Split a TURN-RES Cluster into k Parts (kmeans)}
\description{
	If there is a giant or nuisance component in the clustering (which may have arisen from numerous narrow
	bridges in the dataset), this helper function splits it into k smaller clusters using kmeans.
}

\usage{
clsSplit(x, cls.spec, k)
}
\arguments{
  \item{x}{a \link{cTurn} object containing the clustering that requires splitting.}
  \item{cls.spec}{a vector enumerating the cluster(s) that require splitting.}
  \item{k}{the number of clusters to split apart the specified into.}
}

\details{This is essentially a wrapper for the \link{kmeans} function, but deals with the nuts and bolts of retaining
		the required structure of the input object, \code{x}.}

\value{
	An object of class \code{cTurn}. This will be the same object as \code{x}, but with the \code{$cluster} vector
	now reflecting the desired split.
}

\author{
	Alex Bird, \email{alex.bird@boots.co.uk}
}

\examples{
x = rnorm(30,0,1);
mean(x) #approx 0.5. Nothing to do with this!
}
