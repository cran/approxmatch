\name{approxmatch-package}
\alias{approxmatch-package}
\alias{approxmatch}
\docType{package}
\title{
\packageTitle{approxmatch}.
}
\description{
\packageDescription{approxmatch}
}
\details{

%%The DESCRIPTION file:
%%\packageDESCRIPTION{approxmatch}
\packageIndices{approxmatch}
An R package for creating matched strata with multiple treatments. Default 
design for a stratum structure is one unit from each treatment, but, other 
designs can be specified.  User can also fine match/ near fine match on 
one or more categorical covariates, e.g. sex and age group.

The main functions of the package are \code{kwaymatching} and \code{tripletmatch}.
These functions take as input the distance structure of multiple groups and the grouping
information to create an approximately optimal multigroup design minimizing the 
total distance. A distance structure can be calculated as per requirement by 
the \code{multigrp_dist_struc} function. 

The algorithm used to create matched design is an approximation algorithm developed
by Karmakar, Small and Rosenbaum (2019).  The design built is guaranteed to be 
close to the optimal matched design of the specified structure.

IMPORTANT NOTE: In order to perform matching, \code{kwaymatching} requires the 
user to load the optmatch (>= 0.9-1) package separately. A manual loading is 
required due to software license issues. If the package is not loaded, the 
\code{kwaymatching} command will fail with an error saying the optmatch package 
is not present. Reference to optmatch is given below.
}
\author{
\packageAuthor{approxmatch}

Maintainer: \packageMaintainer{approxmatch}
}
\references{
Crama, Y. and Spieksma, F. C. R. (1992), Approximation algorithms for three-dimensional
assignment problems with triangle inequalities, European Journal of Operational Research
60, 273--279.


Hansen, B.B. and Klopfer, S.O. (2006) Optimal full matching and related designs via network flows,
JCGS 15 609--627.


Karmakar, B., Small, D. S. and Rosenbaum, P. R. (2019) Using Approximation Algorithms to Build 
Evidence Factors and Related Designs for Observational Studies, Journal of Computational 
and Graphical Statistics, 28, 698--709.
}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
%%\seealso{
%%~~ Optional links to other man pages, e.g. ~~
%%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%%}
\examples{
## See kwaymatching for usage
}
