\name{clsTurnRes}
\alias{clsTurnRes}
\alias{cTurn-class}

\title{TURN-RES Clustering}
\description{
	Implementation of the TURN-RES clustering algorithm (Foss, 2002). TURN-RES is a density based clustering 
	algorithm, but achieves superior efficiency and usability over other methods such as DBSCAN. Neighbour
	estimation is achieved through cyclically sorting the dataset over all its dimensions, but note that
	each datapoint is given only one neighbour in every perpendicular direction, and not necessarily the closest one.
}

\usage{
clsTurnRes(data, r, summarise = F, min.size = "Auto", base.cls = "None", 
			phi = 0.8)
}
\arguments{
  \item{data}{\emph{required}. A numeric data frame or matrix where each column is a dimension to be clustered over. 
				Alternatively a cTurn object; ie a previous output of this function.}
  \item{r}{\emph{required}. Resolution parameter for TURN-RES. Think of this like the adjustment wheel of a microscope.
			The smaller the value, the higher the granularity of the clustering. A high resolution will quantize
			data to a coarser grid. The purpose of the \link{clsMRes} function is inform the value of this parameter.}
  \item{summarise}{Output can be summarised if purpose of clustering was for top level metrics. When \link{clsMRes} calls
			\code{clsTurnRes}, it only requires statistics of the run rather than the full output, so a summary is returned.
			There are two levels of summary. (1) is the more highly summarised - simply n, k, %clustered, table of frequencies.
			(2) also returns the cluster vector.}
  \item{min.size}{The minimum size a cluster must be for classification as a cluster. Any clusters smaller than \code{min.size}
			will be considered as noise. The default value is \code{n/100}, so a cluster must represent at least 1\% of the dataset.
			A numeric value must be supplied is interpreted absolutely rather than as a proportion.}
  \item{base.cls}{A somewhat experimental notion. A vector of the known cluster membership of the dataset, such that each row
			of the dataset corresponds to the respective row in the vector, where an integer will specify the cluster number,
			and \code{NA} denotes unknown cluster membership. This is intended to force the algorithm to separate or agglomerate
			clusters based on prior information, but the cutoffs between the forced separation can be an unnatural shape.}
  \item{phi}{Another parameter of the algorithm which determines the density required to agglomerate points. In theory,
			the choice of this parameter is arbitrary (see references), as it effectively 'scales' the resolution parameter.
			There has been no formal proof of this, hence the option to tweak it.}

}

\details{While not completely parameterless, the user is only required to specify \code{r}, the resolution of the clustering. However,
		the algorithm has much more power when paired with its parent function \code{\link{clsMRes}}, which iterates through a sequence
		of values to reveal the structure of the data and aid with parameter selection.

		Desired clusters may be found across different values of the parameter \code{r}. While \code{\link{formaliseClusters}} is designed
		to take in arguments of a \code{\link{clsMR}} object across multiple resolutions, there may be instances where the analyst wants
		to split open a giant cluster for a given resolution. The function \code{\link{clsSplit}} can be called to partition a specified
		cluster(s) into k separate clusters.
		
		In the TURN paper below, a second algorithm , TURN-CUT was proposed to automatically determine the choice of \code{r}. This
		algorithm is in principal similar to the 'elbow method' of determining number of clusters. This has been omitted due to concerns 
		of over-fit and a proposal that an exploratory approach would anyway be preferred. Philosophically, there may be no "best" choice 
		of parameter, as even a given objective may yield a number of different "best" parameters on the same dataset.}

\value{
	An object of class \code{cTurn}. The cluster membership vector can be found in the slot \code{$cluster}. 
	\code{cTurn} objects have a number of generic functions available: \code{print}, \code{summary} and \code{plot}.
}

\note{
	In order to avoid copying the dataset to each cTurn object, instead the name is saved as item \code{$dataset.name}.
	The data will then be retrieved in function calls via \code{get(dataset.name, env = .GlobalEnv)}, which
	means that the user must ensure that the dataset variable name is not changed. This is obviously a suboptimal procedure,
	but given the package is to be used with large datasets, it is also inadvisable to make a copy for every object,
	particularly if dozens of different cluster calls are to be made in quick succession.
}

\references{
\cite{Foss, A. (2002) A Parameterless Method for Efficiently Discover Clusters of Arbitrary Shape in Large Datasets. University of Alberta Canada.}
}

\author{
	Alex Bird, \email{alex.bird@boots.co.uk}
}

\seealso{\code{\link{clsMRes}} for determining the resolution; \code{\link{clsSplit}} for splitting a given cluster into k clusters}

\examples{
#Toy Example
data <- matrix(runif(200),100,2)
cls <- clsTurnRes(data, r = 0.1)

#Cluster Summary
summary(cls)

#Parallel Coordinate Plot
plot(cls)
}

