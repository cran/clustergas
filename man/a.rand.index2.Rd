\name{a.rand.index2}
\alias{a.rand.index2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Adjusted Rand Index.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Adjusted Rand Index measures the extent of agreement between two 
clusterings. When Adjusted Rand Index value increases, the 
agreement grows. 
}
\usage{
a.rand.index2(refpartition, clustering)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{refpartition}{
%%     ~~Describe \code{refpartition} here~~
A reference partition in form of a list of vectors representing 
a clustering. Vector c(0,0) most be added at the end of the list.
}
  \item{clustering}{
%%     ~~Describe \code{clustering} here~~
As the above argument. A list of vectors representing the clustering
that would be compared with the reference partition. Every vector of 
the list represents a cluster of the clustering. Vector c(0,0) most 
be added at the end of the list. If the clustering (as a list) is part 
of a dendrogram, then it is not necessary to add vector c(0,0) at the 
end of the list.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function compares a reference partition of a dataset with a 
clustering using the Adjusted Rand Index as a similarity measure. Note
that the object-data of the dataset are represented as numbers starting 
at 1. Therefore, object-data in each cluster (vector) are represented as
numbers.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A numeric value representing the similarity digree between both 
clusterings.
}
\references{
%% ~put references to the literature/web site here ~
Rand WM: Objective criteria for the evaluation of clustering methods. 
Journal of the American Statistical Association 1971, 66:846-850.
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
\code{\link{degree.similarity}}, \code{\link{rand.index}}, 
\code{\link{jacard.coef}}, \code{\link{minkowski}}
}

\examples{
refpart <- list(c(1, 3, 4), c(5, 2), c(6, 7, 8), c(0, 0))
clust <- list(c(8, 1, 2, 7, 3), c(5, 4, 6), c(0, 0))
a.rand.index2(refpart, clust)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ math }% __ONLY ONE__ keyword per line
