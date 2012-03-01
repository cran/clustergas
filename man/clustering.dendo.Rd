\name{clustering.dendo}
\alias{clustering.dendo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Evolutionary strategy to improve a clustering.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function improves a clustering by exchanging or moving object-data
between the clusters of such a clustering. At the end, it builds a 
dendrogram in the way as crossover operator "crossover3". See function 
\code{\link{agnes.gas}} for more information on this operator.
}
\usage{
clustering.dendo(clustering, level, times = 500, count = 10, 
                 fitness.f = "clustering.fitness", homog.f = "cluster.fitness",
                 print.fitness = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clustering}{
%%     ~~Describe \code{clustering} here~~
Clustering to improve (a list of vectors, each vector represents a cluster).
}
  \item{level}{
%%     ~~Describe \code{level} here~~
Level at which this clustering belongs.
}
  \item{times}{
%%     ~~Describe \code{times} here~~
Exchange number of each iteration. 
}
  \item{count}{
%%     ~~Describe \code{count} here~~
Iteration number.
}
  \item{fitness.f}{
%%     ~~Describe \code{fitness.f} here~~
Name of the fitness function that measures the quality of the dendrograms.
We provide several fitness functions:
\code{\link{fitness.mean}}, \code{\link{sigma.dendo.fitness1}},
\code{\link{fitness.mean.silhouette}}, \code{\link{fitness.meanD.HS}} and
\code{\link{sigma.meanD.HS}}.
}
  \item{homog.f}{
%%     ~~Describe \code{homog.f} here~~
Name of the homogeneity function used by the fitness function of the above
parameter to measure the homogeneity of a cluster. All fitness functions use
homogeneity function "cluster.fitness". Small values of homogeneity
represent clusters of high quality whereas large values represent clusters
of low quality.
}
  \item{print.fitness}{
%%     ~~Describe \code{print.fitness} here~~
If it is TRUE, the fitness values (before and after improving the clustering) 
and the time spent by the function execution are printed on the screen.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list (a dendrogram) of lists (clusterings of the dendrogram) of vectors
(clusters) representing the dendrogram built from the improved clustering. 
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
\note{
%%  ~~further notes~~
Function \code{\link{init.variables}} most be executed before
using any other function of this package. It initializes the
internal variables.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{agnes.gas}}, \code{\link{local.search}},
\code{\link{improve.dendo}}, \code{\link{complete.tree}}, 
\code{\link{dendrogram.graph}}
}
\examples{
library(cluster)
data(votes.repub)
vr.matrix <- as.matrix(daisy(votes.repub, metric = "euclidean", stand = TRUE))
dendro <- agnes.gas(vr.matrix, 5, 10, part = 3/4)
clust <- dendro[[6]] #a clustering of 7 clusters, level 6.
dendro2 <- clustering.dendo(clust, 6, times = 10, count = 20, print.fitness = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ methods }% __ONLY ONE__ keyword per line

