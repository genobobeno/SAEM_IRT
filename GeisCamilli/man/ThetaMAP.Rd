\name{ThetaMAP}
\alias{ThetaMAP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ThetaMAP(aa, bb, cc, rp, settings)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aa}{
%%     ~~Describe \code{aa} here~~
}
  \item{bb}{
%%     ~~Describe \code{bb} here~~
}
  \item{cc}{
%%     ~~Describe \code{cc} here~~
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
function (aa, bb, cc, rp, settings) 
{
    xi <- as.matrix(cbind(aa, bb))
    if (settings$Adim == 1) {
        tgrid <- seq(-3.3, 3.3, length.out = 761)
        p <- ProbOgive(xi = xi, theta = tgrid)
        THAT <- vector()
        for (n in 1:nrow(rp)) {
            y <- t(as.matrix(rp[n, ]) \%*\% t(as.matrix(rep(1, 
                nrow(as.matrix(tgrid))))))
            THAT <- c(THAT, tgrid[which.max(apply(y * log(p) + 
                (1 - y) * log(1 - p), 1, sum))])
        }
    }
    else {
        tgrid <- seq(-3.3, 3.3, length.out = 761)
        THAT <- mat.or.vec(nrow(rp), settings$Adim)
        for (n in 1:nrow(rp)) {
            xy <- rep(1, settings$Adim)
            xy0 <- rep(0, settings$Adim)
            It <- 0
            tsearch <- mat.or.vec(length(tgrid), settings$Adim)
            while (sum(xy != xy0) > 0 | It > 500) {
                xy0 <- xy
                for (q in 1:settings$Adim) {
                  tsearch[, q] <- tgrid
                  p <- ProbOgive(xi = xi, theta = tsearch)
                  y <- t(as.matrix(rp[n, ]) \%*\% t(as.matrix(rep(1, 
                    nrow(as.matrix(tgrid))))))
                  xy[q] <- tgrid[which.max(apply(y * log(p) + 
                    (1 - y) * log(1 - p), 1, sum))]
                  tsearch[, q] <- xy[q]
                  if (!prod(is.finite(xy))) 
                    xy <- xy0 + rep(0.1, settings$Adim)
                }
                It <- It + 1
            }
            THAT[n, ] <- xy
        }
    }
    return(THAT)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
