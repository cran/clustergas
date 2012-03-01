\name{fitness.mean.silhouette}
\alias{fitness.mean.silhouette}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Fitness function measuring the dendrogram silhouette width.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes the dendrogram fitness based on the
mean value of the silhouette width of each clustering 
(\code{\link{silhouette.clust}}) of the one. 
}

\usage{
fitness.mean.silhouette(dendogram)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dendogram}{
%%     ~~Describe \code{dendogram} here~~
A dendrogram represented by a list (dendrogram) of lists
(clusterings of the dendrogram) of vectors (clusters).
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
A value representing the dendrogram fitness.
}
\references{
%% ~put references to the literature/web site here ~
L. Kaufman, P. J. Rousseeuw, Finding Groups in Data.
An Introduction to Clustering Analysis, John Wiley &
Sons, Inc., Hoboken, New Jersey, 2005.
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
\code{\link{fitness.mean}}, \code{\link{sigma.dendo.fitness1}},
\code{\link{fitness.meanD.HS}}, \code{\link{sigma.meanD.HS}}.
}
\examples{
library(cluster)
data(votes.repub)
vr.matrix <- as.matrix(daisy(votes.repub, metric = "euclidean", stand = TRUE))
dendro <- agnes.gas(vr.matrix, 5, 10, part = 3/4)
fitness <- fitness.mean.silhouette(dendro)
fitness
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ optimize }% __ONLY ONE__ keyword per line 

