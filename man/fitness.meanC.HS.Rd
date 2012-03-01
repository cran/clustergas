\name{fitness.meanC.HS}
\alias{fitness.meanC.HS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Clustering fitness function based on separation and homogeneity.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes the clustering fitness based on the
difference between the separation and  homogeneity average 
of the clusters.
}
\usage{
fitness.meanC.HS(clustering)
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
A value representing the clustering fitness.
}
\references{
%% ~put references to the literature/web site here ~
D. Jiang, C. Tang, A. Zhang, Cluster analysis for gene expression data: 
A survey, IEEE Transactions on Knowledge and Data Engineering 16 (11) 
(2004) 1370-1386.
}
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
\code{\link{clustering.fitness}}, \code{\link{sigma.clust.fitness}},
\code{\link{silhouette.clust}}, \code{\link{sigma.clust.HS}}.
}
\examples{
library(cluster)
data(votes.repub)
vr.matrix <- as.matrix(daisy(votes.repub, metric = "euclidean", stand = TRUE))
dendro <- agnes.gas(vr.matrix, 5, 10, part=3/4)
clust <- dendro[[5]]
fitness <- fitness.meanC.HS(clust)
fitness
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ optimize }% __ONLY ONE__ keyword per line