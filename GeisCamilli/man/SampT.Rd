\name{SampT}
\alias{SampT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SampT(aa, bb, zz, rp, prior)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aa}{
%%     ~~Describe \code{aa} here~~
}
  \item{bb}{
%%     ~~Describe \code{bb} here~~
}
  \item{zz}{
%%     ~~Describe \code{zz} here~~
}
  \item{rp}{
%%     ~~Describe \code{rp} here~~
}
  \item{prior}{
%%     ~~Describe \code{prior} here~~
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
function (aa, bb, zz, rp, prior) 
{
    N <- nrow(rp)
    J <- ncol(rp)
    if (ncol(as.matrix(aa)) > 1) {
        A2 <- t(as.matrix(aa)) \%*\% as.matrix(aa)
        V <- ginv(A2 + diag(ncol(as.matrix(aa))))
    }
    else {
        A2 <- as.numeric(t(as.matrix(aa)) \%*\% as.matrix(aa))
        V <- 1/A2
    }
    ZB <- zz + t(as.matrix(bb) \%*\% t(as.matrix(rep(1, N))))
    if (ncol(as.matrix(aa)) > 1) {
        tHat <- t(V \%*\% t(as.matrix(aa)) \%*\% t(ZB))
    }
    else {
        AT <- t(as.matrix(aa) \%*\% t(as.matrix(rep(1, N))))
        tHat <- rowSums((ZB) * AT) * V
    }
    that <- mat.or.vec(N, ncol(as.matrix(aa)))
    if (ncol(as.matrix(aa)) == 1) {
        Mean <- (tHat/V + prior$tmu/prior$tsigma)/(1/V + 1/prior$tsigma)
        SD <- sqrt(1/(1/V + 1/prior$tsigma))
        that <- rnorm(N, mean = Mean, sd = SD)
    }
    else {
        that <- mvrnormArma(N, tHat, V)
    }
    return(that)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
