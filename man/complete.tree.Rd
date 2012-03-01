\name{complete.tree}
\alias{complete.tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Adds the missing levels of a dendrogram to reach its original size.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function complete the levels of a dendrogram whose size was reduced by
parameter "part" (see function \code{\link{agnes.gas}}) to its original size.
}
\usage{
complete.tree(dendogram)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dendogram}{
%%     ~~Describe \code{dendogram} here~~
A dendrogram whose size is less than its original size.
That is, the original size of a dendrogram is the size of
the used dataset less 2.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
To build the down part of the input dendrogram, a divisive strategy that divides 
a cluster (chosen on the last level of the dendrogram) into two ones is executed
in order to form the next lower level. And so on until it reaches the first level
of the new dendrogram. This function is usually used in conjunction with function 
\code{\link{dendrogram.graph}}, which converts a dendrogram to an object of class 
"agnes" (see package "cluster"). Class "agnes" represents an agglomerative 
hierarchical clustering of a dataset.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A dendrogram whose size is the size of the dataset less 2.
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
\code{\link{agnes.gas}}, \code{\link{local.search}},
\code{\link{improve.dendo}}, \code{\link{clustering.dendo}},
\code{\link{dendrogram.graph}}
}
\examples{
library(cluster)
data(votes.repub)
vr.matrix <- as.matrix(daisy(votes.repub, metric = "euclidean", stand = TRUE))
dendro <- agnes.gas(vr.matrix, 5, 10, part=3/4)
length(dendro)
fulldendro <- complete.tree(dendro)
length(fulldendro)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
\keyword{ utilities }% __ONLY ONE__ keyword per line 
