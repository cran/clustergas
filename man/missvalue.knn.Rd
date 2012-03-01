\name{missvalue.knn}
\alias{missvalue.knn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Replacing missing values using k-nearest neighbors.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function replaces missing values by values computed 
from the technique of k-nearest neighbors.
}
\usage{
missvalue.knn(datamat, K = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datamat}{
%%     ~~Describe \code{datamat} here~~
Data matrix (dataset) to remove missing values.
}
  \item{K}{
%%     ~~Describe \code{K} here~~
The number of neighbors used to compute the new
value that replaces the missing value.
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
A data matrix without missing values.
}
\references{
%% ~put references to the literature/web site here ~
Jain A. K., Dubes RC: Algorithms for Clustering Data. 
Prentice Hall Englewood Cliffs, New Jersey 07632 1998.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{standard}}.
}
\examples{
library(cluster)
data(votes.repub)
newdata <- missvalue.knn(votes.repub, K = 5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
\keyword{ utilities }% __ONLY ONE__ keyword per line
