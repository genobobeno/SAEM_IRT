\name{SampZ}
\alias{SampZ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SampZ(aa, bb, that, rp, w = NA, prl = settings$parallel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aa}{
%%     ~~Describe \code{aa} here~~
}
  \item{bb}{
%%     ~~Describe \code{bb} here~~
}
  \item{that}{
%%     ~~Describe \code{that} here~~
}
  \item{rp}{
%%     ~~Describe \code{rp} here~~
}
  \item{w}{
%%     ~~Describe \code{w} here~~
}
  \item{prl}{
%%     ~~Describe \code{prl} here~~
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
function (aa, bb, that, rp, w = NA, prl = settings$parallel) 
{
    Eta <- t(as.matrix(aa) \%*\% t(as.matrix(that))) - t(as.matrix(bb) \%*\% 
        t(as.matrix(rep(1, nrow(rp)))))
    Ez <- mat.or.vec(nrow(rp), ncol(rp))
    if (!is.na(w)) {
        Ez[which(w == 1)] <- TruncNorm(Eta[which(w == 1)], "high", 
            prl)
        Ez[which(w == 0)] <- TruncNorm(Eta[which(w == 0)], "low", 
            prl)
    }
    else {
        Ez[which(rp == 1)] <- TruncNorm(Eta[which(rp == 1)], 
            "high", prl)
        Ez[which(rp == 0)] <- TruncNorm(Eta[which(rp == 0)], 
            "low", prl)
    }
    return(Ez)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
