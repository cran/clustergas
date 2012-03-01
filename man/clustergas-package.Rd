\name{clustergas-package}
\alias{clustergas-package}
\alias{clustergas}
\docType{package}
\title{
A hierarchical clustering method based on genetic algorithms.
}

\description{
This package develops a genetic algorithm to build hierarchical
clusterings from a dataset. It additionally provides a set of 
functions applied to cluster analysis. This package most be 
initialized before its use, to do this, function "init.variables"
most be call to create the variables of internal use. Finally,
the main function of this package is "agnes.gas" (also known as 
HCGA or EMHC) that runs the genetic algorithm to build dendrograms 
from a dataset. A dendrogram of this package is represented as a 
list of lists of vectors. That is, each clustering of the dendrogram 
is represented as list of vectors, where each vector is a cluster 
and each list is a clustering. Hence, all nested clusterings of the 
dendrogram are collected in lists to form a new list representing a 
dendrogram.\cr

Depends: R (>= 2.10.0), cluster
}

\details{
\tabular{ll}{
Package: \tab clustergas\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-03-01\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}
}

\author{
Jose A. Castellanos-Garzon <jantonio_cu@ieee.org>,
Fernando Diaz <fdiaz@infor.uva.es>.\cr

Maintainer: Jose A. Castellanos-Garzon <jantonio_cu@ieee.org>
}

%\references{
%}
\keyword{ package }
%\seealso{
%}
%\examples{
%}
