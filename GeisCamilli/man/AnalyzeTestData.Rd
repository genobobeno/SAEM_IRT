\name{AnalyzeTestData}
\alias{AnalyzeTestData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Begin Analysis of Test Data
}
\description{
This will kick off the process of fitting the data to the Model chosen.
}
\usage{
AnalyzeTestData(RP, settings = settings, verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RP}{
  matrix; This is a response pattern. If dichotomous, 1's and 0's. If polythomous, 0 is still the lowest ordinal category.
  }
  \item{settings}{
  list; This is the list of settings for the algorithm.
  }
  \item{verbose}{
  TRUE/FALSE; This will provide some outputs.
  }
}
\details{
The bulk of this work is based on Stochastic Approximation of the EM Algorithm for Multidimensional Item Response Theory using the Probit link. 
}
\value{
 A list:
\item{ XI }{The generated parameters, slopes and intercept, J X Q + 1 + (guessing==TRUE)+0}
\item{ RP }{The input response pattern}
\item{ THETA}{The generated latent factors}
\item{ A }{ The fitted slopes, J X Q}
\item{ AR}{ The rotated fitted slopes, J X Q}
\item{ B }{The fitted intercepts}
\item{ C }{The fitted guessing parameters}
\item{ xi}{The item-specific fitted matrix of slopes, intercepts, and guessing}
\item{ xiError }{The errors of the parameters in \code{xi}}
\item{ iError}{The Lewis errors of the parameters before SAEM in \code{xi}}
\item{ oError}{The Lewis errors of the parameters during SAEM in \code{xi}}
\item{ gain }{The vector of the gain parameter used at every iteration of the Markov Chain}
\item{ EZ }{The expectation value of the augmented data}
\item{ EZZ }{The expectation of the covariance matrix of the augmented data}
\item{ That }{The EAP of the latent factors, N X Q}
\item{ Tmap }{The MAP of the latent factors, N X Q}
\item{ Tmaprot }{The MAP of the rotated latent factors, N X Q}
\item{ TRmap}{The MAP of the rotated latent factors, N X Q}
\item{ Theta}{The EAP of the rotated latent factors, N X Q}
\item{ Trot=TROT }{}
\item{ EmpSE}{The empirical errors of the parameters in \code{xi}}
\item{ ThetaFix}{The Empirical Bayesian estimate of the latent factors, N X Q}
\item{ settings}{The list of settings for parameters used in the algorithm }
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Eugene Geis
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
