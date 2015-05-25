\name{clsMRes}
\alias{clsMRes}
\alias{clsMR-class}

\title{Perform Clustering for Multiple Resolution Parameters}
\description{
	The initial stage of cluster investigation using nectr. clsMRes shows how the clusters grow 
	and agglomerate across increasing resolutions of the clsTurnRes algorithm. 
}
\usage{
clsMRes(data, r.start = NA, r.max = Inf, keep = FALSE, ...)
}
\arguments{
  \item{data}{\emph{Required}. The data to be clustered. Each column must be a dimension; no labels are permitted.
				\code{data} must be of class \code{data.frame}.}
  \item{keep}{\emph{Optional}. Whether to keep the clustering vector associated with each clsTurnRes iteration.
				By default the algorithm keeps the clustering vector found on each iteration, 
				but one may wish to discard them if the vectors are large.}
  \item{r.start}{\emph{Optional}. For user specified resolution step sequence. This specifies the first number in the sequence.}
  \item{r.max}{\emph{Optional}. For user specified resolution step sequence. This specifies the last number in the sequence - ie the search is terminated if the parameter value is not less than r.max.}
  \item{...}{\emph{Optional} (additional arguments to be passed to clsTurnRes)}
}
\details{Automatic function - should be no need for parameter choice if uncertain. The first goal is to
	find an appropriate sequence of resolutions to run through the clsTurnRes algorithm. A 	starting 
	position is determined automatically, but it may be quicker to specify the starting resolution
	once the analyst is familiar with the dataset.
}

\value{
	An object of class \code{clsMR}. Primary output is to demonstrate structure through increasing resolution
	parameter. Standard visualisation available through the \code{plot} generic function. Object can be 
	interrogated by cluster via the subset [] operator to give details of resolution used.
}

\note{
	In order to avoid copying the dataset to each clsMR object, instead the name is saved as item \code{$dataset.name}.
	The data will then be retrieved in function calls via \code{get(dataset.name, env = .GlobalEnv)}, \emph{which
	means that the user must ensure that the dataset variable name is not changed}. There are obvious problems with this approach
	but the alternative is to waste memory copying potentially large datasets every time the function is called.
}
\author{
	Alex Bird, \email{alex.bird.14@ucl.ac.uk}
}


\seealso{\code{\link{clsTurnRes}}, the clustering algorithm used.}

\examples{
x = rnorm(30,0,1);
mean(x) #approx 0.5. Nothing to do with this!
}
