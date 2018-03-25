\name{WrapX}
\alias{WrapX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Option Propensities Sampler (Polytomous)
}
\description{
Function samples 'eta' option propensities for polytomous test item conditional on current parameters.
}
\usage{
WrapZ(j,A,b,d,theta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{j}{Items for parallelization}
  \item{A}{Current Latent Factor Loadings}
  \item{b}{Current item propensity}
  \item{d}{Previous option propensity}
  \item{theta}{Current Latent Factor}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
matrix of item option propensities
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
