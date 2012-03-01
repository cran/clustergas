\name{clustering.fitness}
\alias{clustering.fitness}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Main fitness function of a clustering.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes the clustering fitness based on the
difference between the separation and  homogeneity mean of the clusters.
}
\usage{
clustering.fitness(clustering)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clustering}{
%%     ~~Describe \code{clustering} here~~
A clustering represented by a list of vectors.
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
A numerical value representing the clustering fitness.
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
\code{\link{sigma.clust.fitness}}, \code{\link{silhouette.clust}}, 
\code{\link{fitness.meanC.HS}}, \code{\link{sigma.clust.HS}}.
}
\examples{
library(cluster)
data(votes.repub)
vr.matrix <- as.matrix(daisy(votes.repub, metric = "euclidean", stand = TRUE))
dendro <- agnes.gas(vr.matrix, 5, 10, part=3/4)
clust <- dendro[[5]]
fitness <- clustering.fitness(clust)
fitness
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ optimize }% __ONLY ONE__ keyword per line
