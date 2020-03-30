\name{covbalance}
\alias{covbalance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check covariate balance of a design.
}
\description{
For a given match, this function evaluates the balance of variables before and
after matching. Balance is evaluated using standardized differences.  
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
covbalance(.data, grouplabel, matches, vars, details)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{.data}{
	A data frame or matrix containing the informations on the \code{vars} for 
	which balance will be checked.
}
  \item{grouplabel}{
	Argument describing the group structure. See the description in the 
	documentation of \code{kwaymatching} function of this package.
}
  \item{matches}{
	A character matrix describing the strata structure of the design. Standard
	use is the output of \code{tripletmatching} or \code{kwaytmatching} function.
	
	Each row of the matrix corresponds to a strata and each entry corresponds to 
	rowname of \code{.data}.
}
  \item{vars}{
	A character vector of the names of the variables for which balance should 
	be checked. 
}
  \item{details}{
	Optional argument. This argument can be used to get other details on the 
	\code{vars} before and after matching. A character vector of names of functions
	which summarizes a vector, e.g. mean. 
}
}
\details{
Standardized difference of the covariates between two groups is computed as 
difference of the means of the variable over the squared root of the average
variance of the variable in the groups. 

For better understanding of the matching, \code{details} can be used. This 
argument can be used to get summaries of the variables before and after matching.
For example, details = c(mean = 'mean', median = 'function(x) quantile(x, probs=.5)') 
given the mean and median of the variables. Only functions that give a single 
number summary can be used! 

Currently, this function cannot be immediately used for a design with different
strata sizes. One way to get around would be to fill in the smaller stratum 
with false units and making all the strata of equal size.
}
\value{
A list consisting of the following elements.
	\item{std_diff}{ Standardized differences of the specified variables before and 
	after matching for every pair of groups.}
	\item{details}{ Only if \code{details} is provided. A list of summaries of the 
	variables before and after matching using the functions specified by \code{details}.}
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Bikram Karmakar 
}
\note{
See \code{\link{kwaymatching}} for usage.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{tripletmatching}}, \code{\link{kwaymatching}}
}
\examples{
	\dontshow{
	data(Dodgeram)
	
	grouplabel = 2*Dodgeram$WITHSABS + 1*Dodgeram$NOSAB + 3*Dodgeram$optSAB
	## An example strata structure
	## Ideally, these should be calculated from the output of the 
	##	\code{kwaymatching} function.  For example, if \code{res}
	##	is the output of this function, use \code{res$matches}.
	
	matches = as.matrix(sample(rownames(Dodgeram), 500), ncol = 5)
	
	vars = c("AGE", "SEX.2", "IMPACT3.3", "DR_DRINK")
	details = c('std_diff')
	
	covbalance(.data=Dodgeram, grouplabel=grouplabel, 
	                 matches = matches, vars = vars, details) }
	data(Dodgeram)
	
	## An example strata structure
	matches = as.matrix(sample(rownames(Dodgeram), 500), ncol = 5)
	
	vars = c("AGE", "SEX.2", "IMPACT3.3", "DR_DRINK")
	details = c('std_diff', 'mean', 'function(x) diff(range(x))', 
	                         'function(x) quantile(x, probs = .9)')
	names(details) <- c('std_diff', 'mean', 'range', '90perc')
	
	covbalance(.data=Dodgeram, grouplabel=c("NOSAB", "optSAB", "WITHSABS"), 
	                 matches = matches, vars = vars, details)
					 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ covariate balance }% use one of  RShowDoc("KEYWORDS")
