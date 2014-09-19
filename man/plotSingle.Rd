\name{plotSingle}
\alias{plotSingle}

\title{Plot Method for a single cluster of a cTurn object}
\description{
	A view of a single cluster within a cTurn object. Default plotting will output all clusters to a panelled Principal
	Component plot. If more detail is required on a particular cluster, this plot wrapper can be used for either a 
	boxplot or a single Principal Component plot.
}

\usage{
plotSingle(x, i, type = "box")
}
\arguments{
  \item{x}{a \link{cTurn} object.}
  \item{i}{the cluster in question. For principal component plots, a vector of clusters can be specified.}
  \item{type}{\emph{optional} The default plot type is a boxplot, but \code{"pcp"} can be specified for a principal
		component plot. If \code{"pcp"} option is used, multiple \code{i} can be specified and viewed on the same plot.
	}
}

\seealso{\link{ggplot2} is used for all plot types.}

\author{
	Alex Bird, \email{alex.bird@boots.co.uk}
}

\examples{
x = rnorm(30,0,1);
mean(x) #approx 0.5. Nothing to do with this!
}
