\name{standard}
\alias{standard}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Normalizing a dataset.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function normalizes a dataset to mean 0 and variance 1.
}
\usage{
standard(mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
%%     ~~Describe \code{mat} here~~
Dataset to normalize.
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
The normalized dataset. 
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{missvalue.knn}}.
}
\examples{
library(cluster)
data(votes.repub)
newdata <- standard(votes.repub)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
\keyword{ utilities }% __ONLY ONE__ keyword per line
