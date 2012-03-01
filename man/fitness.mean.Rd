\name{fitness.mean}
\alias{fitness.mean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Main fitness function of a dendrogram.
}

\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This is the main fitness function of a dendrogram in the package. This 
function computes the dendrogram fitness based on the mean value of 
the fitness of each clustering of the dendrogram. So, the dendrogram 
fitness is computed by the mean fitness of the dendrogram clusterings 
using function \code{\link{clustering.fitness}}.  
}

\usage{
fitness.mean(dendogram)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dendogram}{
%%     ~~Describe \code{dendogram} here~~
A dendrogram represented by a list (dendrogram) of lists 
(clusterings of the dendrogram) of vectors (clusters).
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
A value representing the dendrogram fitness.
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
\code{\link{fitness.meanD.HS}}, \code{\link{sigma.meanD.HS}}.
}

\examples{
library(cluster)
data(votes.repub)
vr.matrix <- as.matrix(daisy(votes.repub, metric = "euclidean", stand = TRUE))
dendro <- agnes.gas(vr.matrix, 5, 10, part = 3/4)
fitness <- fitness.mean(dendro)
fitness
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ optimize }% __ONLY ONE__ keyword per line 
