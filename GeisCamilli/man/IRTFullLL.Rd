\name{IRTFullLL}
\alias{IRTFullLL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
IRTFullLL(xi, rp, settings)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xi}{
%%     ~~Describe \code{xi} here~~
}
  \item{rp}{
%%     ~~Describe \code{rp} here~~
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
function (xi, rp, settings) 
{
    theta <- GetQuad(settings)
    JQ <- mat.or.vec(nrow(theta), 1) + 1
    L <- rep(0, nrow(rp))
    if (tolower(settings$icc) == "ogive") {
        p <- ProbOgive(xi = xi, theta = theta[, 1], settings = settings) - 
            as.matrix(1e-06 * sin(pi/(2 * max(theta[, 1])) * 
                theta[, 1])) \%*\% rep(1, ncol(rp))
    }
    else if (tolower(settings$icc) == "logistic") {
        p <- ProbIRT(xi = xi, theta = theta[, 1], settings = settings)
    }
    for (n in 1:nrow(rp)) {
        X <- t(as.matrix(rp)[n, ] \%*\% t(JQ))
        umL <- X * log(p) + (1 - X) * log(1 - p)
        umL <- exp(rowSums(umL)) * theta[, 2]
        L[n] <- sum(umL)
    }
    print("Full Likelihood done")
    return(L)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
