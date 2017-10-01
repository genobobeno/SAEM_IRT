\name{ProbIRT}
\alias{ProbIRT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ProbIRT(xi, theta, j = NA, Adim = settings$Adim, guess = settings$guess)
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
  \item{Adim}{
%%     ~~Describe \code{Adim} here~~
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
function (xi, theta, j = NA, Adim = settings$Adim, guess = settings$guess) 
{
    if (Adim > 1) 
        print("NOTE: Multidimensional 2PL or 3PL is not implemented")
    if (Adim == 1) {
        if (guess) {
            aa <- as.matrix(xi[, 1])
            bb <- as.matrix(xi[, 2])
            cc <- as.matrix(xi[, 3])
        }
        else {
            aa <- as.matrix(xi[, 1])
            bb <- as.matrix(xi[, 2])
        }
        if (is.na(j)) {
            Z <- outer(theta, as.vector(bb), "-")
            AT <- t(aa \%*\% t(as.matrix(rep(1, length(theta)))))
            if (guess) 
                C <- t(cc \%*\% t(as.matrix(rep(1, nrow(as.matrix(theta))))))
        }
        else {
            AT <- aa[j]
            Z <- (theta - bb[j])
            if (guess) 
                C <- cc[j]
        }
    }
    ifelse(guess, p <- C + (1 - C) * exp(1.7 * AT * Z)/(1 + exp(1.7 * 
        AT * Z)), p <- exp(1.7 * AT * Z)/(1 + exp(1.7 * AT * 
        Z)))
    return(p)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
