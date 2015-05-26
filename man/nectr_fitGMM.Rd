\name{.nectr.fitGMM}
\alias{.nectr.fitGMM}

\title{Fit a (MAP) Estimate of a Gaussian Mixture Model}
\description{
	This function fits a Gaussian Mixture Model (GMM) via the EM algorithm. This is different to most
	GMM implementations in R since it is the MAP estimate over all parameters. The problem of collapsing
	covariances is avoided as well as providing a useful method of generalising clusters discovered with
	TURN-RES (or other prior knowledge). While further efficiencies can be made, this is one of the most
	efficient GMM fitting functions available in R at present. Details of the derivation of the algorithm
	can be found at \url{https://github.com/spoonbill/nectr}. Please note that this is an internal function
	of the \code{nectr} package and is not strictly intended for end user interaction, hence the oblique
	arguments.
}

\usage{
.nectr.fitGMM(params, hyperparams, eps.stop = 0.001, max.iter = 100,
					skip.llh = 1, silent = FALSE)
}
\arguments{
  \item{params}{\emph{required}. An object of class \code{.cTurn.suffStats'} which should contain mixing probabilities, means,
			covariances and references to the data source.}
  \item{centers}{\emph{required}. An object of class \code{.cTurn.suffStats'} which should contain mixing probabilities, means,
			covariances and references to the data source.}
  \item{eps.stop}{Convergence criterion for the log likelihood}
  \item{max.iter}{Maximum number of iterations allowed for the EM steps before the algorithm exits.}
  \item{skip.llh}{Parameter to set when the log likelihood is evaluated. A value of 3 corresponds to evaluating it every 3rd iteration.}
  \item{silent}{By default the function outputs the increase in log likelihood after every step. 
			No output will be shown during EM iterations if \code{silent = TRUE}}
}

\details{Parameters are initialised to the prior values. There is currently no functionality available for random (or other) starts.
		The function is built specifically for use with TURN-RES and random starts are not necessary when the priors are good. EM is
		guaranteed to converge (save for pathologies removed by using MAP estimation). It is not in general known, however, how many
		iterations will be required. One can observe that even differing hyperparameter values can lead to dramatically different
		convergence times. \code{max.iter} should not be necessary for low dimensional datasets - the algorithm is linear in datapoints
		but cubic in dimensionality. Due to the computational overhead of evaluating the log likelihood, this can be skipped, and evaluated
		only one in \code{skip.llh} times. The convergence criterion is linearly interpolated from the evaluations available. For instance,
		if \code{skip.llh = 5}, the previous iteration is deemed to have increased by \code{(llh(0) - llh(-5))/5}, where llh is the log 
		likelihood.\\
		
		The initial parameters can not be easily constructed by the end users, hence the wrappers available for this procedure. It is 
		easy to obtain the parameters from \code{cTurn} by use of the function \code{.nectr.getSuffStats} or one can coerce \code{kMeans}
		objects into \code{cTurn} objects by use of \code{as.cTurn}. In principle, of course, the parameters can be constructed by inspecting
		an existing \code{.cTurn.suffStats'} object or the code mentioned above. Hyperparameters are constructed using \code{.nectr.getHyper}.\\
		
		Infinite values of hyperparameters are permitted in the function. These will simply hold the relevant parameter(s) constant.
		}

\value{
	\item{pi}{The estimated mixing parameters, ie the proportions of the data in each cluster.}
	\item{mu}{The estimated means of the mixture components.}
	\item{S}{The estimated covariance matrices of the mixture components.}
	\item{cluster}{The expected cluster for each datapoint given the component parameters.}
	\item{llh_ascent}{The value of the log likelihood at each iteration (or each \code{skip.llh}th iteration if relevant).}
	\item{iter}{A tuple containing, the number of iterations, and the number of iterations for which the log likelihood was calculated.}
	
	The object is unclassed; there are no further methods available.
}

\note{
	A section of code from \code{\link[mvtnorm:dmvnorm]{mvtnorm}} was used to calculate gaussian density.
}


\author{
	Alex Bird, \email{alex.bird.14@ucl.ac.uk}
}

\seealso{\code{\link{getEMClusters}} for scoring the data into clusters using the parameters estimated from this function.
		\code{\link{formaliseClusters}} for guided use of these functions or if using input from \code{\link{clsMRes}} /
		\code{\link{clsTurnRes}} objects.}

\examples{
x = rnorm(30,0,1);
mean(x) #approx 0.5. Nothing to do with this!
}
