\name{GIFAFullLL}
\alias{GIFAFullLL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GIFAFullLL(aa, bb, zz, tt, prior = prior)
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
  \item{tt}{
%%     ~~Describe \code{tt} here~~
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
function (aa, bb, zz, tt, prior = prior) 
{
    if (length(prior$tmu) < 2) {
        N <- nrow(zz)
        bb <- as.matrix(bb)
        at <- as.matrix(aa) \%*\% matrix(tt, 1, N)
        a2 <- sum(aa^2)
        LL <- (-0.5) * sum(cbind(rowSums(zz^2), as.vector(2 * 
            t(bb) \%*\% t(zz)), sum(bb^2), -2 * colSums(t(zz + 
            as.matrix(rep(1, N)) \%*\% t(bb)) * at), (a2 + 1) * 
            tt^2))
    }
    else {
        N <- nrow(zz)
        bb <- as.matrix(bb)
        at <- as.matrix(aa) \%*\% t(tt)
        LL <- (-0.5) * sum(cbind(rowSums(zz^2), as.vector(2 * 
            t(bb) \%*\% t(zz)), rep(sum(bb^2), N), -2 * colSums(t(zz + 
            as.matrix(rep(1, N)) \%*\% t(bb)) * at), colSums(at^2), 
            rowSums(tt^2)))
    }
    return(LL)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
