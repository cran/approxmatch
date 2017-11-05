\name{Dodgeram}
\alias{Dodgeram}
\docType{data}
\title{
Dodge ram pk 2500 data on side airbag (SAB) usage
}
\description{
This is an example dataset on the description of the Dodge pk 2500 cars involved
in fatal crashes. 
}
\usage{data("Dodgeram")}
\format{
  A data frame with 6953 observations on the following 33 variables.
  \describe{
    \item{\code{X.1}}{a numeric vector.}
    \item{\code{X}}{a numeric vector.}
    \item{\code{indexinSABdata}}{a numeric vector}
    \item{\code{NOSAB}}{a numeric vector. In the earlier period when SAB was not
		available? 1 = yes.}
    \item{\code{WITHoptSAB}}{a numeric vector. Bought SAB in the optional period?
			1 = yes.}
    \item{\code{optSAB}}{a numeric vector. In the optional period when SAB
			was available as an option? 1 = yes.}
    \item{\code{WITHSABS}}{a numeric vector. In the later period when SAB was standard 
			issue? 1 = yes.}
    \item{\code{IMPACT3}}{a numeric vector. Impact code. (0 = no, 1 = right, 2 = behind, 
				3 = left, 4 = front, 9 = other)}
    \item{\code{ROLLOVER1}}{a numeric vector. Rollover occurred? 1 = yes.}
    \item{\code{REST_USE1}}{a numeric vector. Restraint used by the driver?
				1 = yes, 99 = unknown type, 0 = not used.}
    \item{\code{FRpass.REST_USE1}}{a numeric vector. Restraint use by the front 
			right passenger.}
    \item{\code{SP_LIMIT}}{a numeric vector. Speed Limit of the route.}
    \item{\code{AGE}}{a numeric vector. Age of the driver.}
    \item{\code{DR_DRINK}}{a numeric vector. Whether the driver was drinking?
			1 = yes.}
    \item{\code{FR.pass}}{a logical vector. Whether a front right passenger was
			present.}
    \item{\code{FRpass.AGE}}{a numeric vector. Age of the front right passenger
			if present.}
    \item{\code{SEX.2}}{a numeric vector. Sex of the driver. 1 = female.}
    \item{\code{EJECTION.1}}{a numeric vector. Ejection of the driver. 1 = yes.}
    \item{\code{EJECTION.2}}{a numeric vector. Ejection type of driver unknown. 1= yes.}
    \item{\code{IMPACT3.1}}{a numeric vector. Impact from right? 1= yes.}
    \item{\code{IMPACT3.2}}{a numeric vector. Impact from behind? 1= yes.}
    \item{\code{IMPACT3.3}}{a numeric vector. Impact from left? 1= yes.}
    \item{\code{IMPACT3.4}}{a numeric vector. Impact from front? 1= yes.}
    \item{\code{IMPACT3.9}}{a numeric vector. Impact of other type? 1= yes.}
    \item{\code{ROLLOVER1.1}}{a numeric vector. Rollover occurred? 1= yes.}
    \item{\code{FIRE_EXP1.1}}{a numeric vector. Fire occurred in the car? 1= yes.}
    \item{\code{REST_USE1.1}}{a numeric vector. Driver used restraint? 1= yes.}
    \item{\code{REST_USE1.99}}{a numeric vector. Driver restraint use of unknown
			type? 1= yes.}
    \item{\code{FRpass.SEX.2}}{a numeric vector. Sex of the front right passenger.
			2 = female.}
    \item{\code{FRpass.EJECTION.1}}{a numeric vector. Ejection of the front right 
			passenger. 1 = yes.}
    \item{\code{FRpass.EJECTION.2}}{a numeric vector. Ejection type of front right
			passenger unknown. 1= yes.} 
    \item{\code{FRpass.REST_USE1.1}}{a numeric vector. Front right passenger used restraint? 
				1= yes.}
    \item{\code{FRpass.REST_USE1.99}}{a numeric vector. Front right passenger restraint use of unknown
			type? 1= yes.}
  }
}
\details{
Derived from fatal accidents in the USA between 1995 and 2015 recorded by 
Fatality Analysis Reporting System (FARS) of NHTSA. Consists of 
only model years more than 1985. Sample is corrected for selection bias in 
FARS. For further detail of the variables see the code book of FARS.
}
\source{
Fatality Analysis Reporting System of NHTSA, DOT, USA. https://www.nhtsa.gov/research-data/fatality-analysis-reporting-system-fars
}
%%\references{
%%  ~~ possibly secondary sources and usages ~~
%%}
\examples{
data(Dodgeram)
## maybe str(Dodgeram) ; plot(Dodgeram) ...
}
\keyword{datasets}
