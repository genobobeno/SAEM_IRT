\name{NumHessian}
\alias{NumHessian}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
NumHessian(j, xi, rp, pL, settings)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{j}{
%%     ~~Describe \code{j} here~~
}
  \item{xi}{
%%     ~~Describe \code{xi} here~~
}
  \item{rp}{
%%     ~~Describe \code{rp} here~~
}
  \item{pL}{
%%     ~~Describe \code{pL} here~~
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
function (j, xi, rp, pL, settings) 
{
    XI <- xi
    dp <- 0.001
    if (!settings$guess) {
        XI[, 1] <- xi[, 1] - dp
        lowA <- Gradient(j, XI, rp, pL, settings)
        XI[, 1] <- xi[, 1] + dp
        highA <- Gradient(j, XI, rp, pL, settings)
        XI <- xi
        XI[, 2] <- xi[, 2] - dp
        lowB <- Gradient(j, XI, rp, pL, settings)
        XI[, 2] <- xi[, 2] + dp
        highB <- Gradient(j, XI, rp, pL, settings)
        ddA <- 0.5 * (highA - lowA)/dp
        ddB <- 0.5 * (highB - lowB)/dp
        ddL <- rbind(ddA, ddB)
    }
    else {
        XI[, 1] <- xi[, 1] - dp
        lowA <- Gradient(j, XI, rp, pL, settings)
        XI[, 1] <- xi[, 1] + dp
        highA <- Gradient(j, XI, rp, pL, settings)
        XI <- xi
        XI[, 2] <- xi[, 2] - dp
        lowB <- Gradient(j, XI, rp, pL, settings)
        XI[, 2] <- xi[, 2] + dp
        highB <- Gradient(j, XI, rp, pL, settings)
        XI <- xi
        XI[, 3] <- xi[, 3] - dp
        lowC <- Gradient(j, XI, rp, pL, settings)
        XI[, 3] <- xi[, 3] - dp
        highC <- Gradient(j, XI, rp, pL, settings)
        ddA <- 0.5 * (highA - lowA)/dp
        ddB <- 0.5 * (highB - lowB)/dp
        ddC <- 0.5 * (highC - lowC)/dp
        ddL <- rbind(ddA, rbind(ddB, ddC))
    }
    return(ddL)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
