\name{ProbOgive}
\alias{ProbOgive}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ProbOgive(xi, theta, j = NA, guess = settings$guess)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xi}{
%%     ~~Describe \code{xi} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{j}{
%%     ~~Describe \code{j} here~~
}
  \item{guess}{
%%     ~~Describe \code{guess} here~~
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
function (xi, theta, j = NA, guess = settings$guess) 
{
    if (guess) {
        aa <- as.matrix(xi[, 1:(ncol(xi) - 2)])
        bb <- as.matrix(xi[, ncol(xi) - 1])
        cc <- as.matrix(xi[, ncol(xi)])
        ifelse(is.na(j), C <- t(as.matrix(cc) \%*\% t(as.matrix(rep(1, 
            nrow(as.matrix(theta)))))), C <- cc[j])
    }
    else {
        aa <- as.matrix(xi[, 1:(ncol(xi) - 1)])
        bb <- as.matrix(xi[, ncol(xi)])
    }
    if (is.na(j)) {
        AT <- t(aa \%*\% t(as.matrix(theta)))
        Bz <- t(bb \%*\% t(as.matrix(rep(1, nrow(as.matrix(theta))))))
    }
    else {
        Bz <- bb[j]
        ifelse(ncol(aa) > 1, AT <- (aa[j, ] \%*\% t(as.matrix(theta))), 
            AT <- aa[j] * theta)
    }
    if (ncol(aa) == 1) {
        ifelse(!guess, p <- pnorm(AT - Bz), p <- C + (1 - C) * 
            pnorm(AT - Bz))
    }
    else {
        ifelse(!guess, p <- pnorm(AT - Bz), p <- C + (1 - C) * 
            pnorm(AT - Bz))
    }
    return(p)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
