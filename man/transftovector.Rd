\name{transftovector}
\alias{transftovector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Converts a clustering (in the form of list) to a vector structure.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function transforms a clustering represented as 
a list of vectors to a clustering as vector. This last
representation is compatible with the "twins" class 
of the R internal structure. 
}
\usage{
transftovector(clustering)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clustering}{
%%     ~~Describe \code{clustering} here~~
A clustering in the form of a list of vectors.
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
A clustering in the form of a R vector.
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
\code{\link{transf.tree}}.
}
\examples{
library(cluster)
data(votes.repub)
vr.matrix <- as.matrix(daisy(votes.repub, metric = "euclidean", stand = TRUE))
dendro <- agnes.gas(vr.matrix, 5, 10, part = 3/4)
clust <- dendro[[6]] #a clustering of 7 clusters, level 6.
vector.clust <- transftovector(clust)
vector.clust
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
\keyword{ manip }% __ONLY ONE__ keyword per line
