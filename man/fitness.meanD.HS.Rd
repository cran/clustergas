\name{fitness.meanD.HS}
\alias{fitness.meanD.HS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Fitness function of a dendrogram based on 
function "fitness.meanC.HS".
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes the fitness of a dendrogram
based on the difference between homogeneity and 
separation given by function "fitness.meanC.HS".
Which computes the fitness of each clustering of
the dendrogram passed as parameter. Finally, 
"fitness.meanD.HS" computes the mean value of 
all values given by "fitness.meanC.HS".
}

\usage{
fitness.meanD.HS(dendogram)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dendogram}{
%%     ~~Describe \code{dendogram} here~~
A dendrogram represented as a list of lists of vectors.
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
A numerical value representing the fitness value of the dendrogram.
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
\code{\link{sigma.dendo.fitness1}}, \code{\link{fitness.mean.silhouette}},
\code{\link{fitness.mean}}, \code{\link{sigma.meanD.HS}}
}
\examples{
library(cluster)
data(votes.repub)
vr.matrix <- as.matrix(daisy(votes.repub, metric = "euclidean", stand = TRUE))
dendro <- agnes.gas(vr.matrix, 5, 10, part = 3/4)
fitness <- fitness.meanD.HS(dendro)
fitness
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ optimize }% __ONLY ONE__ keyword per line
