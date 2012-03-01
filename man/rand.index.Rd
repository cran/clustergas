\name{rand.index}
\alias{rand.index}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Rand Index measure.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Rand Index measures the extent of agreement between two 
clusterings. When the Rand Index value increases, the 
agreement grows.
}
\usage{
rand.index(nij)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nij}{
%%     ~~Describe \code{nij} here~~
it is a vector returned by function \code{\link{degree.similarity}}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function compares a reference partition of a dataset with a 
clustering using the Rand Index as a similarity measure. Note
that the object-data of the dataset are represented as numbers starting 
at 1. Therefore, object-data in each cluster (vector) are represented as
numbers. The reference partition and the clustering are passed as input
to \code{\link{degree.similarity}}, which returns a vector nij. 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A numeric value representing the similarity digree between two 
clusterings.
}
\references{
%% ~put references to the literature/web site here ~
Halkidi M, Batistakis Y, Vazirgiannis M: On Clustering Validation 
Techniques. In Intelligent Information Systems J. 2001.

Sokal RR: Clustering and Classification: Background and Current 
Directions. Classification and Clustering, Academic Press 1977.
}
%\author{
%%  ~~who you are~~
%}
\note{
%%  ~~further notes~~
\code{\link{init.variables}} function most be executed before using any other
function of this package. It initializes the internal variables.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{degree.similarity}}, \code{\link{a.rand.index2}}, 
\code{\link{jacard.coef}}, \code{\link{minkowski}}
}

\examples{
refpart <- list(c(1, 3, 4), c(5, 2), c(6, 7, 8), c(0, 0))
clust <- list(c(8, 1, 2, 7, 3), c(5, 4, 6), c(0, 0))
nij <- degree.similarity(refpart, clust)
rand.index(nij)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ math }% __ONLY ONE__ keyword per line
