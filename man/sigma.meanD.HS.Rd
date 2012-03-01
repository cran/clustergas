\name{sigma.meanD.HS}
\alias{sigma.meanD.HS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Dendrogram fitness function based on homogeneity and 
separation, checking the cluster size.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This fitness function  is similar to function \code{\link{fitness.meanD.HS}},
but in this case, it adds a new objective (sigma), checking that the number 
of object-data in each cluster of the dendrogram is not very different.
}
\usage{
sigma.meanD.HS(dendogram)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dendogram}{
%%     ~~Describe \code{dendogram} here~~
The dendrogram to compute the fitness.
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
A numerical value representing the dendrogram fitness.
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
\code{\link{fitness.mean}}, \code{\link{fitness.mean.silhouette}},
\code{\link{fitness.meanD.HS}}, \code{\link{sigma.dendo.fitness1}}.
}
\examples{
library(cluster)
data(votes.repub)
vr.matrix <- as.matrix(daisy(votes.repub, metric = "euclidean", stand = TRUE))
dendro <- agnes.gas(vr.matrix, 5, 10, part = 3/4)
fitness <- sigma.meanD.HS(dendro)
fitness
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ optimize }% __ONLY ONE__ keyword per line
