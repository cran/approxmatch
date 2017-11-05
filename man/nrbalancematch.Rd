\name{nrbalancematch}
\alias{nrbalancematch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The is the background function to perform matching using network optimization.
}
\description{
This function takes as input the information on two groups, treated and control,
along with the data structure between the units of these two groups to create a
matching. Custom specification of the categorization of the units for near fine
balanced design and balance constraints can be provided.
}
\usage{
nrbalancematch(cardata.fil, trt_labs, ctrl_labs, stratify, extmatch = NA, 
distmat, balanceValues, nmatch = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cardata.fil}{
	A data frame or matrix of the dataset. The units are identified using the 
	rownames of the data.
}
  \item{trt_labs}{
	A logical vector of same length as nrow(cardata.fil) for the treated units.
}
  \item{ctrl_labs}{
	A logical vector of same length as nrow(cardata.fil) for the control units.
}
  \item{stratify}{
	A character vector of same length as nrow(cardata.fil) specifyting the 
	categorization of the units for fine balancing.
}
  \item{extmatch}{
	The name of the column of cardata.fil on which design will be exactly matched.
	This is optional.
}
  \item{distmat}{
	A numeric matrix of distances of size the number of treated units \eqn{\times}
	number of control units. The rows and columns should be named by the unit 
	names.
}
  \item{balanceValues}{
	A named integer vector of balance values of the categories.
}
  \item{nmatch}{
	Number of control units to be matched to each treated units. A positive 
	integer.
}
}
\details{
Do not use this function directly unless you are certain of the usage. 

IMPORTANT NOTE: In order to perform matching, \code{kwaymatching} requires the 
user to load the optmatch (>= 0.9-1) package separately. The manual loading is 
required due to software license issues. If the package is not loaded the 
\code{nrbalancematch} command will fail with an error saying the optmatch package 
is not present. Reference to optmatch is given below.
}
\value{
A matrix of two columns first column for the treated unit and the second column 
for the control unit. Units in the first column will be repeated \code{nmatch}
times.
}
\references{
Hansen, B.B. and Klopfer, S.O. (2006) Optimal full matching and related designs via network flows,
JCGS 15 609-627.
}
\author{
Bikram Karmakar
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{tripletmatching}}, \code{\link{kwaymatching}}
}
%%\examples{
%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ matching }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
