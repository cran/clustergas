\name{silhouette.clust}
\alias{silhouette.clust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Fitness function measuring the clustering silhouette width.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes the clustering fitness based on the  
silhouette width measure. 
}
\usage{
\method{silhouette}{clust}(\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
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
\code{\link{clustering.fitness}}, \code{\link{sigma.clust.fitness}},
\code{\link{fitness.meanC.HS}}, \code{\link{sigma.clust.HS}}.
}
\examples{
library(cluster)
data(votes.repub)
vr.matrix <- as.matrix(daisy(votes.repub, metric = "euclidean", stand = TRUE))
dendro <- agnes.gas(vr.matrix, 5, 10, part=3/4)
clust <- dendro[[5]]
fitness <- silhouette.clust(clust)
fitness
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ optimize }% __ONLY ONE__ keyword per line
