\name{local.search}
\alias{local.search}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Evolutionary strategy of local search. 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is an evolutionary strategy of local search that improves
the individuals (dendrograms) from a initial population. The initial 
population can randomly be created or read from parameter "init.pop". 
Usually, this function is used to improve the individuals returned 
in a population of individuals by function \code{\link{agnes.gas}}.
}

\usage{
local.search(dist.matrix, pop.size = 10, gen.numb, part = 1/2, 
             mutation = "mutation2", fitness.f = "fitness.mean", 
             homog.f = "cluster.fitness", return.pop = FALSE, 
             init.pop = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dist.matrix}{
%%     ~~Describe \code{dist.matrix} here~~
Proximity matrix (distance matrix) of the dataset to be clustered.
}
  \item{pop.size}{
%%     ~~Describe \code{pop.size} here~~
Size of the population that will be created  or loaded from a file.
}
  \item{gen.numb}{
%%     ~~Describe \code{gen.numb} here~~
Number of generations (iterations) to mate the individuals.
}
  \item{part}{
%%     ~~Describe \code{part} here~~
A fraction representing the part of the dendrograms that
will be removed. By default part = 1/2, which means that the
first half of the dendrograms are removed. In other words,
the size of the dendrograms is reduced to the half. See 
section details of function \code{\link{agnes.gas}} for 
more information. 
}

  \item{mutation}{
%%     ~~Describe \code{mutation} here~~
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
  \item{return.pop}{
If it is FALSE, then this function returns the best individual
of the whole evolutionary process. That is, the best dendrogram.
If it is TRUE, then this function returns a list with two named 
elements. The first element, named "elite", has the best individual
(the best dendrogram) of the whole evolutionary process, and the 
second one, named "lastpop", has the last population of individuals,
generated  in the last generation of the algorithm.
}
\item{init.pop}{
It represents the initial population of individuals (dendrograms)
to be used as a starting point by this function. If "init.pop" is 
NULL, then this function creates a random initial population. 
Otherwise, the initial population of this function is taken from
the argument assigned by the user to "init.pop".
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
This evolutionary strategy mutes all individuals in each generation
by using the passed mutation operator. It only replaces  the 
individuals that improve their fitness. 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
If parameter "return.pop" has FALSE, then this function
returns a list (a dendrogram) of lists (clusterings of the 
dendrogram) of vectors (clusters) representing the best 
dendrogram of the whole evolutionary process. Otherwise, 
it returns a list with two named arguments, where the 
first element is the best dendrogram and second one is 
the last population of individuals as explained for 
parameter "return.pop". Additionally, the time taken by 
the execution of the algorithm and the fitness value of 
the best dendrogram are shown.
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
\code{\link{agnes.gas}}, \code{\link{improve.dendo}}, 
\code{\link{clustering.dendo}}, \code{\link{complete.tree}}, 
\code{\link{dendrogram.graph}}
}

\examples{
library(cluster)
data(votes.repub)
vr.matrix <- as.matrix(daisy(votes.repub, metric = "euclidean", stand = TRUE))

# Example 1
dendro1 <- local.search(vr.matrix, pop.size = 5, 10, part=3/4)

# Example 2
result.list <- agnes.gas(vr.matrix, 5, 10, part = 3/4, return.pop = TRUE)
dendro2 <- result.list$elite
population <- result.list$lastpop
# Improving the individuals of "population".
dendro3 <- local.search(vr.matrix, pop.size = 5, 50, part=3/4, 
                        init.pop = population)
fitness.mean(dendro2)                        
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ methods }% __ONLY ONE__ keyword per line

