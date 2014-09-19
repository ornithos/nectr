\name{clsReNum}
\alias{clsReNum}

\title{Re-number TURN-RES Clustering in line with clsMR Object}
\description{
	Due to the random order of traversal in the TURN-RES algorithm (see \link{clsTurnRes}) and the cumulative nature
	of the \link{clsMRes} function, clusters will be labelled differently between the \link{clsMR} and \link{cTurn} objects.
	In order to facilitate effective inspection of the hierarchy in the \code{clsMR} object, this helper function
	re-numbers the \code{cTurn} object to match the \code{clsMR} object.
}

\usage{
clsReNum(x, y)
}
\arguments{
  \item{x}{a \link{cTurn} object. This is the object to be re-numbered.}
  \item{y}{a \link{clsMR} object. This is the object to extract the desired numeration from.}
}

\details{Naturally the resolution parameter in x must match a resolution in the \code{clsMR} object, otherwise there
		is no like-for-like clustering to match to.}

\value{
	An object of class \code{cTurn}. This will be the same object as \code{x}, but with the \code{$cluster} vector
	re-numbered.
}

\author{
	Alex Bird, \email{alex.bird@boots.co.uk}
}

\examples{
x = rnorm(30,0,1);
mean(x) #approx 0.5. Nothing to do with this!
}
