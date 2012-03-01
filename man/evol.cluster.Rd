\name{evol.cluster}
\alias{evol.cluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Evolutionary strategy to improve a clustering.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is an evolutionary strategy that improves a clustering
according to a fitness function by exchanging two objects between 
two different clusters or moving one object from a cluster to another. 
}
\usage{
evol.cluster(clustering, gen.numb, fitness.f = "clustering.fitness", 
             mop = "gene.change1", print.fitness = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clustering}{
%%     ~~Describe \code{clustering} here~~
The clustering to improve.
}
  \item{gen.numb}{
%%     ~~Describe \code{gen.numb} here~~
Number of generations. Number of times that the operator 
is applied.
}
  \item{fitness.f}{
%%     ~~Describe \code{fitness.f} here~~
Fitness function of the clustering.
}
  \item{mop}{
%%     ~~Describe \code{mop} here~~
Type of operator to use, that can be
"gene.change1" or "gene.change2".
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
Operator "gene.change1" moves one object from a cluster to 
another cluster. Operator "gene.change2" exchanges two objects
located in different clusters.
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
\code{\link{exec.strag}}, \code{\link{clustering.dendo}},
\code{\link{improve.dendo}}.
}
\examples{
library(cluster)
data(votes.repub)
vr.matrix <- as.matrix(daisy(votes.repub, metric = "euclidean", stand = TRUE))
dendro <- agnes.gas(vr.matrix, 5, 10, part = 3/4)
clust <- dendro[[6]] #a clustering of 7 clusters, level 6.
clust2 <- evol.cluster(clust, 100, mop = "gene.change2", print.fitness = TRUE) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ methods }% __ONLY ONE__ keyword per line
