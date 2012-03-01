\name{silhouette.mean}
\alias{silhouette.mean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Clustering silhouette measure. 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function measures silhouette with of a clustering,
relating homogeneity and separation.
}
\usage{
\method{silhouette}{mean}(\dots)
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
A numerical value representing the silhouette 
width of the clustering.
}
\references{
%% ~put references to the literature/web site here ~
Kaufman L, Rousseeuw PJ: Finding Groups in Data. 
An Introduction to Clustering Analysis. John Wiley &
Sons, Inc., Hoboken, New Jersey 2005.
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
\code{\link{H.ave}}, \code{\link{S.ave}}.
}
\examples{
library(cluster)
data(votes.repub)
vr.matrix <- as.matrix(daisy(votes.repub, metric = "euclidean", stand = TRUE))
dendro <- agnes.gas(vr.matrix, 5, 10, part = 3/4)
clust <- dendro[[6]] #a clustering of 7 clusters, level 6.
silhouette.mean(clust)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ math }% __ONLY ONE__ keyword per line
