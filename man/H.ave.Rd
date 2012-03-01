\name{H.ave}
\alias{H.ave}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Clustering homogeneity measure.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function measures the homogeneity of a clustering.
}
\usage{
H.ave(clustering)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clustering}{
%%     ~~Describe \code{clustering} here~~
A clustering in form of a list of vectors.
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
A numerical value representing the homogeneity of
the clustering.
}
\references{
%% ~put references to the literature/web site here ~
Jiang D, Tang C, Zhang A, Cluster Analysis for Gene 
Expression Data: A Survey. IEEE Transactions on 
Knowledge and Data Engineering 2004, 16(11):1370-1386.
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
\code{\link{S.ave}}, \code{\link{silhouette.mean}}.
}
\examples{
library(cluster)
data(votes.repub)
vr.matrix <- as.matrix(daisy(votes.repub, metric = "euclidean", stand = TRUE))
dendro <- agnes.gas(vr.matrix, 5, 10, part = 3/4)
clust <- dendro[[6]] 
H.ave(clust)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ math }% __ONLY ONE__ keyword per line
