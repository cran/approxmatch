\name{kwaymatching}
\alias{kwaymatching}
\alias{tripletmatching}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create approximately optimal matched strata of multiple groups.
}
\description{
This function takes as input a distance structure and grouping labels of units 
to create an strata structure of the units that minimizes the sum distance.  It 
also provides features for near fine balance on one or more nominal 
variables and exact match on a nominal variable.
}
\usage{
kwaymatching(distmat, grouplabel, design, indexgroup = 1, .data, finebalanceVars, 
                       exactmatchon, ordering, reorder = FALSE, verbose = TRUE)
		
tripletmatching(distmat, grouplabel, design, indexgroup = 1, .data, finebalanceVars, 
                       exactmatchon)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distmat}{
	An output of \code{multigrp_dist_struc}. This a named list of matrices. Each 
	list element correspond to the distance matrix between two groups. If group
	labels are 1:k. Then the elements are named '1-2', '1-3',\eqn{\cdots}, 'i-j', \eqn{\cdots}.
	For matching of three groups k=3. If groups are labeled by some character
	string then list elements are named accordingly. For example, if the groups
	are 'grpA', 'grpB', 'grpC'; the elements are 'grpA-grpB', 'grpA-grpC', 
	'grpB-grpB'. 
	
	Each element of the list is a numeric matrix of distances between units of the
	corresponding groups with number of rows equal to the size of the first group,
	and number of columns the size of the second group. The units are identified
	by their names in row or column names of the list elements.
}
  \item{grouplabel}{
	This argument is used to provide information about the group structure of the
	units. There are a few options on how this information can be provided.
	
	1) A numeric vector or a categorical vector. By providing a numeric vector 
	which labels each the units to corresponding groups.
	
	2) A matrix or data frame of dummy variables with number of groups 
	as the number of rows (one for each group) and number of units as the 
	number of columns.
	
	3) A character vector of variable names. If \code{.data} is provided 
	then only name of the variable which contains the information about the 
	grouping can be provided. If the grouping information is encoded in dummy 
	variable, one can also just provide the names of the dummy variables.

	If either the first or the second kind of information is provided, it is 
	expected that the unit are identifiable from the names (in the first case)
	or rownames (in the second case) of \code{grouplabel}.
}
  \item{design}{
	A vector of positive integers specifying the design. If not provided then 
	the default is one unit from each group.
}
  \item{indexgroup}{
	The number or the name of the group to be considered as the index group. 
	The design size of the index group should be 1. Further, for each 
	non-index group the size of the group should be at least the product of 
	design size of that group and size of the index group.
}
  \item{.data}{
	Optional argument but recommended. The data frame or matrix of the dataset.
	The units are recognized by rownames. This is used when further design 
	structure of near fine balance and exact match on certain variable is 
	imposed.
}
  \item{finebalanceVars}{
	An optional character vector of names of the columns of \code{.data} on which the 
	matching will be near fine balanced.
}
  \item{exactmatchon}{
	Name of the column of \code{.data} on which matching will be exact. (optional)
}
  \item{ordering}{
	Optional vector of size the number of groups, specifying the order in which 
	groups will be matched sequentially.
}
  \item{reorder}{
	Optional logical argument for whether matching will be done after permuting 
	the groups randomly.
}
  \item{verbose}{
	A logical argument. If true some details about the implementation may be 
	prompted when running.
}
}
\details{
Only required arguments are distmat, grouplabel, design and indexgroup. \code{.data}
is suggested but not required, in which case units are assumed to be labeled 
1:length(grouplabel). 

Argument \code{.data} must be provided if structure of fine balance and/or exact 
match is imposed by arguments \code{finebalanceVars} and/or \code{exactmatchon}.


IMPORTANT NOTE: In order to perform matching, \code{kwaymatching} requires the 
user to load the optmatch (>= 0.9-1) package separately. The manual loading is 
required due to software license issues. If the package is not loaded the 
\code{kwaymatching} command will fail with an error saying the optmatch package 
is not present. Reference to optmatch is given below.
}
\value{
	A list consisting of the following two elements.
	\item{matches}{ A character matrix of size, size of the indexgroup\eqn{\times}
	\code{sum(design)}. Each row corresponds to a strata and cell values are 
	the units in the strata. }
	\item{cost}{ The cost of the final matching calculated as the sum of 
	the average distances between each pair of groups within every strata, and
	summing over all strata.}
}
\references{
Hansen, B.B. and Klopfer, S.O. (2006) Optimal full matching and related designs via network flows,
JCGS 15 609--627.
}
\author{
Bikram Karmakar
}
\note{
For theoretical guarantee it is expected that the distance structure satisfies 
triangle inequality. If triangle inequality is satisfied than the cost of 
the solution of this algorithm is at most twice that of the optimal match.
Note, for more than two groups the problem of finding the optimal solution is NP-hard.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{	
	\dontshow{
	data(Dodgeram) # dodge ram pk 2500
	Dodgeram <- Dodgeram[sample(rownames(Dodgeram), 1000),]
	grouplabel = 2*Dodgeram$WITHSABS + 1*Dodgeram$NOSAB + 3*Dodgeram$optSAB
	
	# distance components consists of log propensity 
	# distance and rank based Mahalanobis distance.
	components <- list(prop = c("AGE", "IMPACT3"), mahal = c("SEX.2", "AGE", "FIRE_EXP1.1"))
	wgts <- c(10, 5)
	distmat <- multigrp_dist_struc(Dodgeram, grouplabel = grouplabel, components, wgts)	}
	
	## USAGE 1
	\dontrun{
	library(optmatch)
	## User is required to install and load the optmatch package separately,
	
	data(Dodgeram) # dodge ram pk 2500
	
	grouplabel = 2*Dodgeram$WITHSABS + 1*Dodgeram$NOSAB + 3*Dodgeram$optSAB
	
	# distance components consists of log propensity 
	# distance and rank based Mahalanobis distance.
	components <- list(prop = c("AGE", "IMPACT3"), mahal = c("SEX.2", "AGE", "FIRE_EXP1.1"))
	wgts <- c(10, 5)
	distmat <- multigrp_dist_struc(Dodgeram, grouplabel = grouplabel, components, wgts)	
	
	# Matching 
	design = c(1,1,3) # 3 units from the optional period, 1 each from other periods
	indexgroup = 2
	res = tripletmatching(distmat = distmat, grouplabel = grouplabel, design = design, 
	                            indexgroup = indexgroup)
	
	# covariance balance
	details = 'mean'
	covbalance(Dodgeram, grouplabel=c("NOSAB", "optSAB", "WITHSABS"), matches = res, 
	                          vars = c("AGE", "SEX.2", "IMPACT3.3", "DR_DRINK"), details)
	}
	
	## USAGE 2
	\dontrun{
	library(optmatch)
	## User is required to install and load the optmatch package separately,
	
	data(Dodgeram)
	
	# Example distance structure
	components <- list(prop = c("AGE", "SEX.2", "FR.pass", "REST_USE1", "ROLLOVER1",
	            "IMPACT3", "SP_LIMIT", "DR_DRINK", "FIRE_EXP1.1"), 
	            mahal = c("SEX.2", "AGE", "SP_LIMIT", "DR_DRINK"), 				
	            mahal = c("IMPACT3", "REST_USE1"))
	wgts <- c(5, 8, 20)
	
	distmat <- multigrp_dist_struc(Dodgeram, grouplabel = c("NOSAB", "optSAB", "WITHSABS"), 
	                                   components, wgts)
	
	# Matching with fine balance and exact match
	indexgroup = "WITHSABS"
	finebalanceVars = c("ROLLOVER1.1", "FIRE_EXP1.1")
	exactmatchon = "FR.pass"
	
	res = tripletmatching(distmat = distmat, grouplabel = c("NOSAB", "optSAB", "WITHSABS"), 
	                design = c(3,3,1), indexgroup = indexgroup, .data = Dodgeram, 
	                finebalanceVars = finebalanceVars, exactmatchon = exactmatchon)
	
	# covariance balance
	vars = c("AGE", "SEX.2", "IMPACT3.3", "DR_DRINK")
	details = c('std_diff', 'mean', 'function(x) diff(range(x))', 
	                         'function(x) quantile(x, probs = .9)')
	names(details) <- c('std_diff', 'mean', 'range', '90perc')
	
	covbalance(.data=Dodgeram, grouplabel=c("NOSAB", "optSAB", "WITHSABS"), 
	                 matches = res, vars = vars, details)
	
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ matching }% use one of  RShowDoc("KEYWORDS")
\keyword{ design }% __ONLY ONE__ keyword per line
