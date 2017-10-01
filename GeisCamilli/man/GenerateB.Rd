\name{GenerateB}
\alias{GenerateB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GenerateB(j, bdist, bparams)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{j}{
%%     ~~Describe \code{j} here~~
}
  \item{bdist}{
%%     ~~Describe \code{bdist} here~~
}
  \item{bparams}{
%%     ~~Describe \code{bparams} here~~
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
function (j, bdist, bparams) 
{
    if (bdist == "norm") {
        b <- rnorm(j, bparams[1], bparams[2])
    }
    else if (bdist == "unif") {
        b <- runif(j, bparams[1], bparams[2])
    }
    return(b)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
