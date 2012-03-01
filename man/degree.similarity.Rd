\name{degree.similarity}
\alias{degree.similarity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Builds the vector of agreement and disagreement between 
the two clusterings passed as arguments. 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function builds a vector of four components, where
each one has the number of object pairs together (or not) 
in both clusterings.
}
\usage{
degree.similarity(clustering, groundtruth)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clustering}{
%%     ~~Describe \code{clustering} here~~
The clustering to be compared with the other argument.
}
  \item{groundtruth}{
%%     ~~Describe \code{groundtruth} here~~
A clustering which acts as a reference partition of 
the used data set.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The vector returned by this function is used to 
measure the agreement and disagreement through
functions "rand.index", "jacard.coef" and "minkowski".
Therefore, this vector is the input to these functions.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A vector representing in each component, the number of object 
pairs together or separated into each cluster of both clusterings.
}
\references{
%% ~put references to the literature/web site here ~
Jiang D, Tang C, Zhang A: Cluster Analysis for Gene Expression Data: 
A Survey. IEEE Transactions on Knowledge and Data Engineering 2004, 
16(11):1370-1386.
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
\code{\link{rand.index}}, \code{\link{jacard.coef}}, 
\code{\link{minkowski}}.
}
\examples{
library(cluster)
data(votes.repub)
vr.matrix <- as.matrix(daisy(votes.repub, metric = "euclidean", stand = TRUE))
dendro <- agnes.gas(vr.matrix, 5, 100, part=3/4)
clust <- dendro[[6]]
dendro <- agnes.gas(vr.matrix, 5, 100, part=3/4)
ref.part <- dendro[[7]]
nij <- degree.similarity(clust, ref.part)
nij
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }
\keyword{ cluster }% __ONLY ONE__ keyword per line
