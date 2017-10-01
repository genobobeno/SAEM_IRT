\name{InitializeParams}
\alias{InitializeParams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
InitializeParams(rp, settings = settings)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
function (rp, settings = settings) 
{
    J <- ncol(rp)
    N <- nrow(rp)
    if (settings$initialize == "random") {
        if (settings$Adim == 1) {
            THat <- rnorm(N, 0, 1)
            A <- runif(J, 0.5, 1.5)
        }
        else {
            A <- mat.or.vec(J, settings$Adim)
            THat <- mat.or.vec(N, settings$Adim)
            for (i in 1:settings$Adim) {
                A[, i] <- runif(J, 0.5, 1.5)
                THat[, i] <- rnorm(N, 0, 1)
            }
        }
        B <- rnorm(J, 0, 1)
        XI <- cbind(A, B)
        if (settings$guess) {
            C <- runif(J, 0.01, 0.3)
            XI <- cbind(XI, C)
        }
    }
    else if (settings$initialize == "best") {
        if (settings$Adim == 1) {
            A <- rep(1, J)
            THat <- qnorm(0.005 + 0.99 * (rowSums(rp)/J))
        }
        else {
            if (settings$Adim == 1) {
                A0 = 1
            }
            else {
                A0 = matrix(c(rep(0.85, J), rep(0.45, J * (settings$Adim - 
                  1))), J, settings$Adim) + 0.3 * matrix(runif(J * 
                  settings$Adim), J, settings$Adim)
            }
            A <- mat.or.vec(J, settings$Adim) + A0
            THat <- as.matrix(qnorm(0.005 + 0.99 * (rowSums(rp)/J)))
            for (i in 2:settings$Adim) {
                THat <- cbind(THat, THat[, 1])
            }
        }
        B <- qnorm(1 - colSums(rp)/N) * 2
        XI <- cbind(A, B)
        if (settings$guess) {
            C <- runif(J, 0.01, 0.3)
            XI <- cbind(XI, C)
        }
    }
    return(list(XI = XI, THat = THat))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
