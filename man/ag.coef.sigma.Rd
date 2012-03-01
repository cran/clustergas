\name{ag.coef.sigma}
\alias{ag.coef.sigma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Agglomerative coefficient for function "sigma.clust.fitness".
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes the level whose clustering has the maximum
fitness value of the dendrogram according to fitness function
"sigma.clust.fitness".
}
\usage{
ag.coef.sigma(dendogram)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dendogram}{
%%     ~~Describe \code{dendogram} here~~
A dendrogram given by a list of lists of vectors.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function applies "sigma.clust.fitness" to each
clustering of the dendrogram to obtain the level where
the fitness is maximum. This function spends less runtime
than function "ag.coef1".
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A number representing the level with the highest fitness 
of the dendrogram.
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
\code{\link{ag.coef1}}, \code{\link{ag.coef.mean1}}
}
\examples{
library(cluster)
data(votes.repub)
vr.matrix <- as.matrix(daisy(votes.repub, metric = "euclidean", stand = TRUE))
dendro <- agnes.gas(vr.matrix, 5, 10, part = 3/4, fitness.f = "sigma.dendo.fitness1")
ag.coef.sigma(dendro)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ methods }% __ONLY ONE__ keyword per line

