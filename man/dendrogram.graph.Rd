\name{dendrogram.graph}
\alias{dendrogram.graph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Converts a full dendrogram to an object of class "agnes".
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function converts a full dendrogram (with all levels) to an 
object of class "agnes". 
}
\usage{
dendrogram.graph(dendogram, datam, fitness.clust = "clustering.fitness", 
                 homog.f = "cluster.fitness")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dendogram}{
%%     ~~Describe \code{dendogram} here~~
 The dendrogram with all levels to convert it to an object
 of class "agnes". 
}
  \item{datam}{
%%     ~~Describe \code{datam} here~~
used dataset.
}
  \item{fitness.clust}{
%%     ~~Describe \code{fitness.clust} here~~
Name of the fitness function that measures the quality of a clustering.
We provide several fitness functions:
\code{\link{clustering.fitness}}, \code{\link{sigma.clust.fitness}},
\code{\link{silhouette.clust}}, \code{\link{fitness.meanC.HS}}, 
\code{\link{sigma.clust.HS}}.
}
  \item{homog.f}{
%%     ~~Describe \code{homog.f} here~~
Name of the homogeneity function used by the fitness function of the above
parameter to measure the homogeneity of a cluster. All fitness functions use
homogeneity function "cluster.fitness". Small values of homogeneity
represent clusters of high quality whereas large values represent clusters
of low quality.
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
The dendrogram is converted to an object of class "agnes", which
allows to use R functions that handle this type of object. The objects of class 
"agnes" represent an agglomerative hierarchical clustering of a dataset.
Namely, we can use function "plot" that displays the graph of the 
dendrogram on a dataset. See package "cluster" of R for more information on 
"agnes.object".
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
An object of class "agnes" representing a dendrogram.
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
\code{\link{improve.dendo}}, \code{\link{clustering.dendo}},
\code{\link{complete.tree}}, \code{\link{transf.tree}}.
}
\examples{
library(cluster)
data(votes.repub)
vr.matrix <- as.matrix(daisy(votes.repub, metric = "euclidean", stand = TRUE))
dendro <- agnes.gas(vr.matrix, 5, 10, part=3/4)
agn <- dendrogram.graph(complete.tree(dendro), votes.repub)
plot(agn)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs }
\keyword{ cluster }% __ONLY ONE__ keyword per line

