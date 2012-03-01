\name{transf.tree}
\alias{transf.tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Transforms a dendrogram of class "twins" to a dendrogram
as a list of lists of vectors.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function transform a dendrogram of class "twins" to 
a list of lists of vectors. This last dendrogram structure is 
the one used by function "agnes.gas" and all functions of this 
package.
}
\usage{
transf.tree(agn, part = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{agn}{
%%     ~~Describe \code{agn} here~~
A dendrogram of the structure of class "twins".
}
  \item{part}{
%%     ~~Describe \code{part} here~~
If we want to reduce the length of the returned dendrogram,
"part" states the proportion of the dendrogram that is 
removed. For example, "part = 4/5" means that the forth-fifth
of the levels in the dendrogram (from the first level) 
will be removed. See this parameter in function "agnes.gas" 
for more information.
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
A dendrogram as a list of lists of vectors,
internal structure used in this package.
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
\code{\link{transftovector}}, \code{\link{dendrogram.graph}}.
}
\examples{
 library(cluster)
 data(votes.repub)
 agn1 <- agnes(votes.repub, metric = "manhattan", stand = TRUE)
 dendro <- transf.tree(agn1, part = 3/4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
\keyword{ manip }% __ONLY ONE__ keyword per line
