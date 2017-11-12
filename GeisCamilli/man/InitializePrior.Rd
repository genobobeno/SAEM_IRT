\name{InitializePrior}
\alias{InitializePrior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Initialize Values for Model Fit
}
\description{
This function randomizes or makes a best guess for initial values for the model fit. For example, the starting value of the first dimension of each examinee's latent factor is calculated from their overall test score if we use "best" initialization. If we choose "random", the latent factor is sampled from the prior. 
}
\usage{
InitializePrior(rp, settings = settings)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rp}{
matrix; the response pattern, rows = number of examinees and columns = items.
}
  \item{settings}{
list; the configuration parameters for the fit.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list;
\item{XI}{The matrix of initialized fit Parameters; cbind(A,b) if a 2PNO, cbind(A,b,c) if a 3PNO}
\item{THat}{The matrix of initialized latent factors of examinees}
%% ...
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
Structure = CheckParams()
Gen.Data = GenerateTestData(j=10,n=500,structure=Structure)
InitializePrior(Gen.Data$RP,settings = CheckParams(generate=FALSE))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
