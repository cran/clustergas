\name{agnes.gas}
\alias{agnes.gas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Genetic algorithm for data hierarchical clustering.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This is the main function of the package. This function runs a
genetic algorithm that is a hierarchical clustering method building
dendrograms from a dataset. It builds or reads from parameter
"init.pop", an initial population of individuals (dendrogram 
population) to mate them a given number of generations  
through of the used genetic operators.
}

\usage{
agnes.gas(dist.matrix, pop.size = 10, gen.numb = 10, part = 1/2,
          fitness.f = "fitness.mean", homog.f = "cluster.fitness",
          crossover = "crossover2", mutation = "mutation2",
          p.cross = 0.4, p.mut = 0.1, return.pop = FALSE, 
          init.pop = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dist.matrix}{
%%     ~~Describe \code{dist.matrix} here~~
Proximity matrix (distance matrix) of the dataset to be clustered.
}
  \item{pop.size}{
%%     ~~Describe \code{pop.zise} here~~
Size of the population that will be created  or loaded from a file.
}
  \item{gen.numb}{
%%     ~~Describe \code{gen.numb} here~~
Number of generations (iterations) to mate the individuals.
}
  \item{part}{
%%     ~~Describe \code{part} here~~
A fraction representing the part of the dendrograms that
will be removed. By default "part = 1/2", which means that the
first half of the dendrograms is removed. In other words,
the size of the dendrograms is reduced to the half.
}
  \item{fitness.f}{
%%     ~~Describe \code{fitness.f} here~~
Name of the fitness function that measures the quality of the dendrograms of
the genetic algorithm. We provide several fitness functions:
\code{\link{fitness.mean}}, \code{\link{sigma.dendo.fitness1}},
\code{\link{fitness.mean.silhouette}}, \code{\link{fitness.meanD.HS}} and
\code{\link{sigma.meanD.HS}}.
}
  \item{homog.f}{
%%     ~~Describe \code{homog.f} here~
Name of the homogeneity function used by the fitness function of the above
parameter to measure the homogeneity of a cluster. All fitness functions use
homogeneity function "cluster.fitness". Small values of homogeneity
represent clusters of high quality whereas large values represent clusters
of low quality.
}
  \item{crossover}{
%%     ~~Describe \code{crossover} here~~
Name of the crossover operator used by the genetic algorithm. We provide
three types of operators, "crossover", "crosssover2" and "crossover3",
all of them cross two dendrograms to produce a new dendrogram.

The "crossover" operator chooses randomly a level on two dendrograms
and the best clusters (according to homogeneity) of each dendrogram in
the chosen level are taken out to form a new clustering. Which is used
to build the child dendrogram. An agglomerative strategy that randomly merges
two clusters of the current clustering in order to form the next upper levels is used
to build the top part of the dendrogram. A divisive strategy that divides a
cluster (randomly chosen) of the current level into two clusters to form the
next lower levels is used to build the down part of the dendrogram.

The "crossover2" operator performs as the "crossover" one, but in this case,
to build the child dendrogram from the new clustering, the agglomerative strategy
repeatedly merges two clusters (usually, the 30\% of all combinations of two clusters
taken of the clustering are selected) and the most homogeny merge of them is taken out
to form the next upper levels. On the other hand, the divisive strategy that builds
the down part of the dendrogram just divides the less homogeneous cluster to form
the next lower level, and so on.

The "crossover3" operator performs as the "crossover2" one, but in this case, before
building the child dendrogram from the new clustering, an improvement strategy  of
this clustering which exchanges object-data between the clusters of this one is
executed. Usually, this evolutionary strategy carries out 200 exchanges between
randomly selected clusters.
}
  \item{mutation}{
%%     ~~Describe \code{mutation} here~~
Name of the mutation operator used by the genetic algorithm. We provide
three types of operators, "mutation", "mutation2" and "mutation3",
all of them mute a dendrogram to produce another dendrogram.

The "mutation" operator chooses randomly a level of the dendrogram, afterwards,
it merges two randomly selected clusters of the current level to form the
next upper level and so on, until it reaches the top of the new dendrogram.

The "mutation2" operator performs as the "mutation" one. In this case,
it repeatedly merges two randomly selected clusters (usually, the 30\% of all
combinations of two clusters are selected) of the chosen level. The most
homogenous merge is chosen to form the next upper level, and so on, until it
reaches the top of the new dendrogram.

The "mutation3" operator exchanges two object-data between the two clusters of
the penultimate level of the dendrogram (the level that has only two clusters). 
Afterwards, that change is propagated to the lower levels of the dendrogram to 
obtain the new one.
}
  \item{p.cross}{
%%     ~~Describe \code{p.cross} here~~
Crossover likelihood applied to the crossover operator.
}
  \item{p.mut}{
%%     ~~Describe \code{p.mut} here~~
Mutation likelihood applied to the mutation operator.
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
Argument "part" is used to define the new size of the dendrograms in
the whole evolutionary process of the method. Thus, it reduces the
level number of the dendrograms to a size defined by the user. That is, 
"dendrogram size = n - 2 - |part * n|", where
n is the size of the used dataset. That is, if the size of a
dataset is 20 objects, then a dendrogram of this one will have
20 levels. But if for example, we want to delete the fourth part of
the dendrograms, that is, part = 1/4, then 
"20 - 2 - |1/4 * 20| = 13" is the new size of the
dendrograms.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
If parameter "return.pop" has FALSE, then this function returns a list 
(a dendrogram) of lists (clusterings in the dendrogram) of 
vectors (clusters) representing to the best dendrogram of the 
whole evolutionary process in the genetic algorithm. Otherwise, 
it returns a list with two named arguments, where the first 
element is the best dendrogram and second one is the last 
population of individuals as explained for parameter "return.pop".
Additionally, the time taken by the execution of  the algorithm 
and the fitness value of the best dendrogram are shown.
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
\note{
%%  ~~further notes~~
Keep in mind that returning the last population of invidividuals
through parameter "return.pop = TRUE", allows to the user 
running this function (the genetic algorithm) over again but
taking as initial population the last population returned by 
the above running, through parameter "init.pop". This way, the
current solutions can be improved. Moreover, for that purpose,
the last population of individuals generated by "agnes.gas" can 
also be given as initial population to function 
\code{\link{local.search}}, which performs an evolutionary local 
search. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{local.search}}, \code{\link{improve.dendo}}, 
\code{\link{clustering.dendo}}, 
\code{\link{complete.tree}}, \code{\link{dendrogram.graph}}
}
\examples{
library(cluster)
data(votes.repub)
vr.matrix <- as.matrix(daisy(votes.repub, metric = "euclidean", stand = TRUE))
# Example 1
dendro <- agnes.gas(vr.matrix, 5, 10, part = 3/4)

# Example 2
result <- agnes.gas(vr.matrix, 5, 10, part = 3/4, return.pop = TRUE)
best.ind <- result$elite
dendro1 <- agnes.gas(vr.matrix, 5, 50, part = 3/4, init.pop = result$lastpop)
fitness.mean(best.ind)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ methods }% __ONLY ONE__ keyword per line
