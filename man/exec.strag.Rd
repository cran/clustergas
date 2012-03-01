\name{exec.strag}
\alias{exec.strag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Evolutionary strategy to improve a clustering.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is an evolutionary strategy that improves a clustering
according to a fitness function by applying two genetic operators that
exchanging two objects between two different clusters and moving one 
object from a cluster to another.
}
\usage{
exec.strag(clustering, times = 1000, count = 10, fitness.f = "clustering.fitness", 
           print.fitness = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clustering}{
%%     ~~Describe \code{clustering} here~~
The clustering to improve.
}
  \item{times}{
%%     ~~Describe \code{times} here~~
Number of times that each operator is applied.
}
  \item{count}{
%%     ~~Describe \code{count} here~~
Number of times that the "times" parameter is applied.
That is, number of times that both operators are applied
each "times" times.
}
  \item{fitness.f}{
%%     ~~Describe \code{fitness.f} here~~
Fitness function applied to the clustering.
}
  \item{print.fitness}{
%%     ~~Describe \code{print.fitness} here~~
If it is TRUE, the fitness values of the initial
and final clustering  as well as the runtime are
screen printed.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The operators used are "gene.change1" and "gene.change2", which
move one object from a cluster to another and exchange two objects
located in different clusters, respectively.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A new clustering as a result of the applied evolution process.
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
\code{\link{evol.cluster}}, \code{\link{clustering.dendo}},
\code{\link{improve.dendo}}.
}
\examples{
library(cluster)
data(votes.repub)
vr.matrix <- as.matrix(daisy(votes.repub, metric = "euclidean", stand = TRUE))
dendro <- agnes.gas(vr.matrix, 5, 10, part = 3/4)
clust <- dendro[[6]] #a clustering of 7 clusters, level 6.
clust2 <- exec.strag(clust, times = 100, count = 5, print.fitness = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ methods }% __ONLY ONE__ keyword per line
