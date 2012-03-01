\name{init.variables}
\alias{init.variables}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
 Initialization of the global variables of the package.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function creates and initializes all global variables used by
the functions defined in this package. It most be call
before using any other function.
}
\usage{
init.variables(dist.matrix, pop.size = 10, gen.numb = 10, p.cross = 0.4, 
               p.mut = 0.1, part = 1/2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dist.matrix}{
%%     ~~Describe \code{dist.matrix} here~~
Proximity matrix of the dataset to be clustered.
}
  \item{pop.size}{
%%     ~~Describe \code{pop.zise} here~~
Size of the population to be used by the genetic algorithm functions.
}
  \item{gen.numb}{
%%     ~~Describe \code{gen.numb} here~~
Number of generations (iterations) of the genetic algorithm.
}
  \item{p.cross}{
%%     ~~Describe \code{p.cross} here~~
Value of the crossover likelihood of the genetic algorithm.
}
  \item{p.mut}{
%%     ~~Describe \code{p.mut} here~~
Value of the  mutation likelihood of the genetic algorithm.
}
  \item{part}{
%%     ~~Describe \code{part} here~~
It is a fraction representing the part of the dendrogram
that we want to remove (example: part = 3/4). 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Argument "part" is used to remove the firts part of the 
dendrogram if we do not want that it has all its
levels. That is, if a dendrogram has 20 levels and part= 2/3,
then it really will have 20 - 2 - |2/3 * 20| = 5 levels, 
removing the 15 first levels. By default "part = 1/2", which 
means that the first half of the dendrogram is removed.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
No return value.
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
\note{
%%  ~~further notes~~
This is the firts function that most be called before using any other
function different of "agnes.gas". Because it initializes 
important global variables.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{agnes.gas}} and \code{\link{local.search}}.
}
\examples{
library(cluster)
data(votes.repub)
vr.matrix <- as.matrix(daisy(votes.repub, metric = "euclidean", stand = TRUE))
init.variables(vr.matrix, part = 3/4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
\keyword{ IO }% __ONLY ONE__ keyword per line
