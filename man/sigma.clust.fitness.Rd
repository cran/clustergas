\name{sigma.clust.fitness}
\alias{sigma.clust.fitness}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Clustering fitness function considering the size of the clusters.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This fitness function  is similar to function \code{\link{clustering.fitness}},
but in this case, it adds a new objective (sigma) that measures the number of 
data-object of each clusters. 
}
\usage{
sigma.clust.fitness(clustering)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clustering}{
%%     ~~Describe \code{clustering} here~~
A clustering represented by a list of vectors.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function additionally checks that the size of the clusters
is not very different from each other. 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A value representing the clustering fitness.
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
\code{\link{clustering.fitness}}, \code{\link{silhouette.clust}},
\code{\link{fitness.meanC.HS}}, \code{\link{sigma.clust.HS}}.
}
\examples{
library(cluster)
data(votes.repub)
vr.matrix <- as.matrix(daisy(votes.repub, metric = "euclidean", stand = TRUE))
dendro <- agnes.gas(vr.matrix, 5, 10, part=3/4)
clust <- dendro[[5]]
fitness <- sigma.clust.fitness(clust)
fitness
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ optimize }% __ONLY ONE__ keyword per line
