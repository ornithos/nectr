## nectr - Exploratory Clustering

A number of years ago I spent a while researching density-based clustering algorithms for exploratory analysis. I was looking to find regions of high density within customer attributes, which might inform a wider campaign or strategy. k-means wasn't really appropriate since it isn't directly concerned with concentration of density, instead returning spherical/voronoi clusters which are optimally spread out. Indeed under a uniform distribution, k-means will happily return an arbitrary partition with no indication of its inappropriateness.

Parametric mixture models (most famously Gaussian Mixture Models) are more flexible with cluster shape, and confer various advantages against k-means, but they still assume too much of the generating distribution. My experience of commercial data has not shown much evidence of any parametric distribution that I am familiar with. Various non-parametric alternatives exist such as DBSCAN or mean-shift, but these are highly sensitive to parameter specification, and suffer from quadratic time or worse. Another approach is Spectral Clustering, using the graph Laplacian, but unless approximated is at least quadratic, and one must specify the number of clusters in advance. 

The algorithm I chose is one called TURN-RES from a [paper](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.7.1966&rep=rep1&type=pdf) by Foss \& Zaïane in 2002. It is an approximation of a density search -- more from a compsci perspective than a statistical one -- which runs efficiently (in `O(nlogn)` time) and having essentially no parameters. It approximates density estimation by quantizing the space into a discrete grid, and computing the distance to all neighbours in axis-oriented directions. It is therefore less effective for highly correlated distributions. The implementation of the main algorithm -- packaged in R -- is described below, and follows the paper except for a few unclear details for which I have taken my best guess.

### TURN-RES
Two complementary algorithms are proposed in the Foss paper, TurnCut and TURN-RES. I believe TURN is an acronym, but I can't help you out with its meaning. The idea is that TURN-RES is the workhorse clustering algorithm, and TurnCut is its partner for choosing a good parameter. The latter proceeds in sequential manner calling TURN-RES and seeking to find the 'elbow' or 'knee' of various metrics (similar to those used for k-means). The authors justify the 'parameter free' description by use of this second algorithm. However, while the results for the experiments in the paper appear to be pretty good for TurnCut, I have not found it particularly useful in certain real world datasets. This is in part due to the subjectivity inherent in clustering; often a continuum of plausible results can be found for a particular parameter. Also, one often finds that different clusters have different natural resolutions and a global resolution parameter is doomed to failure. 

This observation motivated a more hands-on approach, whereby instead of an algorithm attempting to find a solution to an ill-posed problem, the user is permitted to view a summary of intermediate results. By inspection, a suitable value or values of the resolution can be chosen. This process also can be a useful tool for the analyst, in that understanding the density at various resolutions allows one to build up a more complete picture of the underlying distribution.

### Algorithmic overview

We'll quickly explore the algorithm via a little animation. A window of some 2D data is shown below. The original data is shown in grey.

![dead](/_images/clsTurnRes.gif "steps of TURN-RES")

+ Impose a grid at the given resolution (this is the only parameter choice in the algorithm).
+ Quantize the dataset -- move each datapoint to the nearest point on the grid.
+ For each datapoint, consider the neighbours in all axis-aligned directions.
+ If sufficient neighbours exist, the datapoint is labelled an 'interior' point (shown here in red) and may be connected to each of its neighbours.
+ The union of all neighbouring interior points and their (exterior) neighbours form a cluster.
+ The data is returned to the original non-quantised space, together with the cluster labels.

The algorithm can be written to execute in `O(dnlogn)` time, by partitioning the data into rows, sorting by column position, and obtaining horizontal neighbours. Vertical neighbours are discovered analogously. For higher dimensions, this process continues sequentially as sorts across each dimension. Thus datapoints are neighbours if they are distance 1 away *in at least 1 coordinate direction* with respect to the current quantisation constant. My experience is using only 5-6 dimensions; I haven't explored its performance for medium-large dimensional space.

## Parsimonious model compression using the GMM
While the non-parametric ideal is to be admired, there are many attributes that are less desirable. Assigning unseen, test-time data to a cluster is a challenge. It is also by definition difficult to interpret a cluster, since they may be arbitrary shape and span a wide range of different positions in the sample space. For the customer segmentation case we also want to assign outliers with very sparse density to their closest cluster, although it is useful to know that they are unusual. In order to help in these situations, a Gaussian Mixture Model has also been implemented in the package. It is integrated sufficiently well that the user need not know virtually anything about their construction. Because the exploratory work has already been performed, we have chosen the number of clusters K, and the TURN clusters may be used as priors<sup>1</sup>. By selecting the strength of the prior, one can interpolate between the extremes of a mixture of Gaussians with the moments of the TURN clusters, or simply an initialisation of the GMM from these cluster centers. Derivations of the MAP estimates are available in the repo.

<br>
--------------------------
# nectr: Non-Parametric Exploratory Clustering using TURN-RES

These functions are bundled together in the R-package `nectr` and can be found [here](https://github.com/spoonbill/nectr). The following is a simple demonstration of the algorithm on some dummy data. First, let's create some:

## Create a 5 cluster synthetic dataset in 3 dimensions:

```R
library(MASS)
n <- c(175, 60, 50, 135, 105) * 1000
noise.n <- 75000

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

fakeData <- as.data.frame(rbind(mvrnorm(n[1], mu[[1]], sigma[[1]]), mvrnorm(n[2], mu[[2]], sigma[[2]]),
                            mvrnorm(n[3], mu[[3]], sigma[[3]]), mvrnorm(n[4], mu[[4]], sigma[[4]]), 
                            mvrnorm(n[5],mu[[5]], sigma[[5]]), matrix(runif(noise.n*3,-10,10), noise.n , 3)))
names(fakeData) <- c("X", "Y", "Z")
```

We have also added some uniformly distributed noise to make the problem slightly more challenging. The dataset is shown below - we have n = 600,000 datapoints (not all shown!), and d = 3 dimensions. The colour plot shows the ground truth, the black-and-white plot is the data as received by the algorithm.

![fake data](/_images/clusters_together.png "fake data")


## Tuning the resolution parameter:
The core algorithm TURN-RES requires one parameter - the quantization resolution. The lower the resolution, the more neighbours a datapoint has, which increases the cluster sizes. As discussed above, searching for the optimal resolution is not automated, but the `clsMRes` function scans through some likely values in advance for you:

```R
library(nectr)
multires <- clsMRes(fakeData, keep = TRUE)
```
![alt text](/_images/Rmultires.PNG "R output")

The clsMRes function iteratively reduces the resolution until no clusters are visible, and then increases until most datapoints are put in the same cluster. The `keep = TRUE` argument tells the function to keep all of the cluster assignments generated by the various iterations of TURN-RES. We may then use these to build a hierachical picture of how the clusters agglomerate. These 8 calls to the TURN-RES algorithm collectively took 30 seconds on a 2.7GHz i5 (single) processor, and in principle could be trivially parallelised.


## Exploring the dataset via a hierarchical representation
Inspecting the agglomeration schedule, that is, how the clusters agglomerate as the scale is increased:

```R
plot(multires)
```

![alt text](/_images/agglom.png "plot clsMR object")

This plot is key to exploring the dataset with `nectr`. The clusters are plotted in order of increasing resolution on the x axis. In this example, 5 clusters were found in the first run of TURN with resolution r = 0.102, numbered 1:5. This corresponds to the lowest level resolution, and demonstrates that there are only 5 highly concentrated areas of density. (The properties of which can be explored, as shown in the the following section). Five (larger) clusters were discovered with the next resolution of r = 0.144, and the agglomeration schedule shows the clusters which are merged into a superset of themselves - for example cluster 6 is a superset of cluster 1. We can infer that cluster 6 contains a region of high density, attached or surrounded by a region of slightly lower density. At the next iteration, decreasing the resolution again, clusters 6, 9 and 10 are merged into cluster 12. At this stage we may wish to view clusters 6, 9 and 10 to check whether it makes sense to merge these clusters. And so the process continues. The below animation shows the clusters discovered over increasing resolution.

![alt text](/_images/clusterAnim.gif "Clusters over increasing resolution")

## Investigating discovered clusters further

```R
plot(multires, c(8,12,13))
```
![alt text](/_images/pcp8_12_13.png "Parallel Coordinate Plot")

The code specifies that we wish to see samples from clusters 8, 12 and 13 as indexed by the cluster agglomeration plot above. Sufficient samples are usually plotted to have an implicit view of the mean. The plots use [Parallel Coordinate plots](https://en.wikipedia.org/wiki/Parallel_coordinates), where succesive dimensions are plotted on the x-axis. Note that clusters do not have to be of the same resolution to be plotted together, although it often makes better sense to do so. In the above plot, the order of the graphs corresponds to the order in which the clusters were specified. In particular, cluster 12 (plot 2) demonstrates 3 quite different signatures, suggesting that merging clusters 6, 9 and 10 may not be desired. In the below we plot all the clusters at the previous level motivated by the undesirable composition of cluster 12.

```R
plot(multires,6:10)
```
![alt text](/_images/pcp_6_10.png "Parallel Coordinate Plot")

Here we see 5 distinct clusters, and indeed the plots have very similar means to those specified by our dataset generation code. At this stage we may decide that 5 is a good choice for K (the number of clusters), and for each cluster choose the highest level that it is not agglomerated with another. Further calls to TURN-RES with user-defined resolutions may be necessary, which can be called via: 
```R
modelTurnRes <- clsTurnRes(data, r= ...)
```
where the resolution is specified by `r=...`. Good values may be estimated given those already used in the printout from `clsMres(...)`. Sometimes we may wish to stop here.


### Creating the summary GMM model
Now if a parametric form of this exploratory phase is required, we can pass the output to a Gaussian Mixture Model, centered at the clusters discovered using TURN. The function `clsSpecifyModel` takes as arguments the multiresolution object created by `clsMRes(...)`, the cluster numbers which we wish to keep in the GMM model (again these do not have to be from the same resolution level), and a guess at the percentage of noise. The output of `clsMres(...)` may help with this. The GMM will use the specified clusters as priors, and if a `noise.pct` is specified, an additional large variance noise cluster is admitted into the GMM to prevent outliers skewing the results<sup>2</sup>.

```R
gmmSpec <- clsSpecifyModel(multires, clusters = c(6,7,8,9,10), noise.pct = 0.1)
```
Now we have specified the model, we must fit it. The following function takes the model as an argument, and `alpha`, `beta`, which are the strength of the priors used. These are scaled such that each prior strength should be specified in the interval `[0,1]`. `0` allows the GMM fit complete freedom in fitting, and the resulting model may look very different from the TURN model, at the other extreme, using `1` will simply return the covariance and means of the existing TURN clusters. Any value in the middle is a trade off between these objectives, and below we have used `0.2`. Different values may be specified for each cluster if desired.

```R
modGMM <- clsGMM(gmmSpec, alpha = 0.2, beta = 0.2)
```
The model is fitted by standard EM, which in this example took about 25 seconds for me. The results are shown below; both the 5 fitted clusters, discovered using TURN, and described using the GMM, and in the second, the effect of filtering out low density observations as noise.

![final result](/_images/clustersGMM_together.png "final result")

<br>
##### Footnotes:
<sup>1</sup> the GMM uses the MAP not MLE estimate. This also avoids the collapsing covariance issue common to MLE.

<sup>2</sup> Gaussian distributions are notoriously sensitive to outliers, due to the double exponential decay of the density.

---------------------------------------------------------------------------------------
### Original README notes for R package

This is an R package intended for exploratory analysis of large N, small D datasets such as are typically used by industry analysts. It is based around an implementation of the non-parametric clustering method authored by Foss, 2002 (TURN-RES). The benefits are that no structure is assumed a priori and that the algorithm asymptotically scales linearithmically in N and D. Because we must perform a search on the tuning parameter, we can build a picture of how clusters merge, or are agglomerated as the resolution is decreased. This can then be used to inform a final more parsimonious clustering or dimensionality reducing technique.


#### Notes on TURN-RES implementation:
*  TURN-RES algorithm works fast and effectively. It may be slightly different to the paper mentioned above, due to a couple of details which have required some guesswork. 
*  Algorithm complexity is O(nlogn) dominated by 2*D calls to a mergesort implementation (D = dimension).
*  The algorithm suffers more than many from the curse of dimensionality. Because search is effectively performed conditioned on a constant value of all dimensions bar the given one, the number of datapoints for a set of (D-1) given parameters falls exponentially in D. By way of example, I would not recommend using a dataset of size less than 500,000 for D greater than 5.
* Quantizing - the critical step required to reduce nearest neighbour search to an ordering problem - does induce unintuitive cluster boundaries at times, particularly for a low resolution quantization. 
* Bridges between clusters of similar density can easily be formed through noise. Unlike DBSCAN, we cannot tune by requiring a greater number of points to form the bridge. Therefore excluding toy datasets, a 'best' parameter cannot be found or even defined. Discovering some clusters will naturally be balanced by unintentionally merging others.  Similar issues arise when different variance scale clusters exist in the same dataset
* While the original paper showed good results with the dedicated TURN-CUT algorithm to chose the resolution parameter, I have found real world data to be less accomodating, and due to poor results have omitted this tuning stage from the implementation.
  + All errors or omissions from (Foss, 2002) or indeed anywhere else remain my own.

Other algorithms implemented:
------------------------------------------
TURN-RES is a useful algorithm for data exploration but is fairly useless as a dimensionality reduction method. The Gaussian Mixture Model is included in order to house a more effective model for this purpose in which we can use the results of TURN-RES as a prior.
* Gaussian Mixture Model fitted with EM (MAP estimate used to avoid the collapsing variance issue).

