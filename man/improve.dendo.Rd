\name{improve.dendo}
\alias{improve.dendo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Evolutionary strategy to improve a dendrogram.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function carries out an improvement of a dendrogram by 
using a mutation operator.
}
\usage{
improve.dendo(dendo, ite = 10, operator = "mutation2", 
              fitness.f = "fitness.mean", homog.f = "cluster.fitness", 
              print.fitness = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dendo}{
%%     ~~Describe \code{dendo} here~~
Dendrogram to improve (as a list of lists of vectors).
}
  \item{ite}{
%%     ~~Describe \code{ite} here~~
Iteration number of the strategy.
}
  \item{operator}{
%%     ~~Describe \code{operator} here~~
Name of the mutation operator used by the evolutionary strategy.
We provide three types of operators, "mutation", "mutation2"
and "mutation3", all of them mute a dendrogram to produce another
dendrogram. These operators are explained in function \code{\link{agnes.gas}}.
}
  \item{fitness.f}{
%%     ~~Describe \code{fitness.f} here~~
Name of the fitness function that measures the quality of the dendrograms. 
We provide several fitness functions:
\code{\link{fitness.mean}}, \code{\link{sigma.dendo.fitness1}},
\code{\link{fitness.mean.silhouette}}, \code{\link{fitness.meanD.HS}} and
\code{\link{sigma.meanD.HS}}.
}
  \item{homog.f}{
%%     ~~Describe \code{homog.f} here~~
Name of the homogeneity function used by the fitness function of the above
parameter to measure the homogeneity of a cluster. All fitness functions use
homogeneity function "cluster.fitness". Small values of homogeneity
represent clusters of high quality whereas large values represent clusters
of low quality.
}
  \item{print.fitness}{
%%     ~~Describe \code{print.fitness} here~~
If it is TRUE, the fitness values (before and after improving the input dendrogram) 
and the time spent by the function execution are printed on the screen.
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
A list (a dendrogram) of lists (clusterings of the dendrogram) of vectors
(clusters) representing a dendrogram as a improvement of the input dendrogram. 
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
\code{\link{clustering.dendo}}, \code{\link{complete.tree}}, 
\code{\link{dendrogram.graph}}
}
\examples{
library(cluster)
data(votes.repub)
vr.matrix <- as.matrix(daisy(votes.repub, metric = "euclidean", stand = TRUE))
dendro <- agnes.gas(vr.matrix, 5, 10, part = 3/4)
dendro2 <- improve.dendo(dendro, ite = 20, print.fitness = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ methods }% __ONLY ONE__ keyword per line

