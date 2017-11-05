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

The main functions of the package are \code{kwaymatching} and \code{tripletmatch}.
These functions take as input the distance structure of multiple groups and the grouping
information to create an approximately optimal multigroup design minimizing the 
total distance. A distance structure can be calculated as per requirement by 
the \code{multigrp_dist_struc} function. 

IMPORTANT NOTE: In order to perform matching, \code{kwaymatching} requires the 
user to load the optmatch (>= 0.9-1) package separately. The manual loading is 
required due to software license issues. If the package is not loaded the 
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
