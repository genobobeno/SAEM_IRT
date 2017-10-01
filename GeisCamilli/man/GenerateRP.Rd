\name{GenerateRP}
\alias{GenerateRP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GenerateRP(xi, theta, structure = structure)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xi}{
%%     ~~Describe \code{xi} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{structure}{
%%     ~~Describe \code{structure} here~~
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
function (xi, theta, structure = structure) 
{
    RP <- mat.or.vec(nrow(as.matrix(theta)), nrow(as.matrix(xi)))
    if (tolower(structure$icc) == "ogive") {
        P <- ProbOgive(xi, theta, guess = structure$guess)
    }
    else if (tolower(structure$icc) == "logistic") {
        P <- ProbIRT(xi, theta, guess = structure$guess)
    }
    for (i in 1:nrow(RP)) for (j in 1:ncol(RP)) RP[i, j] <- rbinom(1, 
        size = 1, prob = P[i, j])
    return(RP)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
