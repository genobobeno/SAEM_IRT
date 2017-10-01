\name{AnalyzeTestData}
\alias{AnalyzeTestData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
AnalyzeTestData(RP, settings = settings)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RP}{
%%     ~~Describe \code{RP} here~~
}
  \item{settings}{
%%     ~~Describe \code{settings} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
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

## The function is currently defined as
function (RP, settings = settings) 
{
    stopwatch <- Sys.time()
    stopifnot(nrow(RP) > ncol(RP))
    rp <- RP
    J <- ncol(rp)
    N <- nrow(rp)
    print(paste("*****************   Analyzing Test Data, J =", 
        J, ";  N =", N, "; Output =", settings$estfile, "  ******************"))
    Init <- InitializeParams(rp, settings = settings)
    if (tolower(settings$model) == "gifa") {
        Estimates <- GoGIFA(rp, init = Init, settings = settings)
    }
    else if (tolower(settings$model) == "irt") {
        Estimates <- GoIRT(rp, init = Init, settings = settings)
    }
    else {
        print(paste("Model", settings$model, "not implemented yet"))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
