\name{multigrp_dist_struc}
\alias{multigrp_dist_struc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct the distance structure for the multiple groups.
}
\description{
This function can be used to calculate the distance structure for multiple 
groups. The output of this function can be feed into the argument \code{distmat}
of the main functions \code{kwaymatching} and \code{tripletmatching}.
}
\usage{
multigrp_dist_struc(.data, grouplabel, components, wgts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{.data}{
	The data frame or matrix of the dataset. 
}
  \item{grouplabel}{
	The information on the group structure of the units. See description of 
	\code{\link{kwaymatching}} for details on the argument.
}
  \item{components}{
	A list specifying the components of the distance structure.
	Each element of the list is a character vector of column names
	of the .data on which 'distance' will be calculated. The 
	element names specify the function to be used to calculate to 
	distance of two groups. Element named
	'prop' indicates the propensity distance where the propensity
	is calculated from the specified variable. 
	Element named 'mahal' or 'Mahalanobis' for rank based
	Mahalanobis distance. 
	
	User can spacify their own distance function. For example, a function 
	\code{myDist} should be a function of two arguments:
	a logical vector of the first group indicator and a data matrix. It
	should return a numeric matrix of size number of units of 
	first group\eqn{\times}number of units of 
	second group. See details for an example.
	
}
  \item{wgts}{
	A non-negative numeric vector of weights of the components. 
}
}
\details{
This function can be used to get distance structure suitable 
for creating the distances between the units of the groups. 

For an example of the kind of user defined distance function that can be used see \code{smahal} below.

}
\value{
A list describing the distance structure. For detail see the description of 
the argument \code{distmat} in the function \code{\link{kwaymatching}}.
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Bikram Karmakar
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{kwaymatching}}, \code{\link{tripletmatching}}
}
\examples{
	data(Dodgeram)
	\dontshow{	Dodgeram <- Dodgeram[sample(rownames(Dodgeram), 1000),] }
	# Example distance structure
	components <- list(prop = c("AGE", "SEX.2", "FR.pass", "REST_USE1", "ROLLOVER1",
	            "IMPACT3", "SP_LIMIT", "DR_DRINK", "FIRE_EXP1.1"), 
	            mahal = c("SEX.2", "AGE", "SP_LIMIT", "DR_DRINK"), 				
	            mahal = c("IMPACT3", "REST_USE1"))
	wgts <- c(5, 8, 20)
	
	distmat <- multigrp_dist_struc(Dodgeram, 
            grouplabel = c("NOSAB","optSAB","WITHSABS"), components, wgts)
	
	## Propensity score caliper can be implemented mannually
	
	distmat <- multigrp_dist_struc(Dodgeram, 
            grouplabel = c("NOSAB","optSAB","WITHSABS"), 
            list(mahal = c("SEX.2", "AGE", "SP_LIMIT", "DR_DRINK"), 				
                mahal = c("IMPACT3", "REST_USE1")), wgts=c(2, 5))
	distmat_prop <- multigrp_dist_struc(Dodgeram, 
                grouplabel = c("NOSAB", "optSAB", "WITHSABS"), 
              list(prop = c("AGE", "SEX.2", "FR.pass", "REST_USE1", "ROLLOVER1",
	            "IMPACT3", "SP_LIMIT", "DR_DRINK", "FIRE_EXP1.1")), 1)
		
	## Distance structure with caliper
	for(i in 1:length(distmat))
			distmat[[i]][distmat_prop[[i]]>.2] <- 100*max(distmat[[i]])
			
				
									   
## An example function for argument detail.

smahal <-  function(z,X){
        X<-as.matrix(X)
		n<-dim(X)[1]
		rownames(X)<-1:n
		k<-dim(X)[2]
		m<-sum(z)
		for (j in 1:k) X[,j]<-rank(X[,j])
		cv<-cov(X)
		vuntied<-var(1:n)
		rat<-sqrt(vuntied/diag(cv))
		cv<-diag(rat)%*%cv%*%diag(rat)
		out<-matrix(NA,m,n-m)
		Xc<-X[z==0,,drop=FALSE]
		Xt<-X[z==1,,drop=FALSE]
		rownames(out)<-rownames(X)[z==1]
		colnames(out)<-rownames(X)[z==0]
		#library(MASS)
		icov<-ginv(cv)
		for (i in 1:m) 
			out[i,]<-mahalanobis(Xc,Xt[i,],icov,inverted=T)
		sqrt(out)
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distance structure }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
