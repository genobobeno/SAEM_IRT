\name{Gradient}
\alias{Gradient}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Gradient(j, xi, rp, pL, settings)
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
    theta <- GetQuad(settings)
    NQ <- mat.or.vec(nrow(rp), 1) + 1
    JQ <- mat.or.vec(nrow(theta), 1) + 1
    if (settings$guess) {
        aa <- xi[, 1]
        bb <- xi[, 2]
        cc <- xi[, 3]
    }
    else {
        aa <- xi[, 1]
        bb <- xi[, 2]
    }
    if (tolower(settings$icc) == "ogive") {
        Z <- NQ \%*\% t(aa[j] * theta[, 1] - bb[j])
        TT <- NQ \%*\% t(theta[, 1])
        p <- pnorm(aa[j] * theta[, 1] - bb[j]) - 1e-06 * sin(pi/(2 * 
            max(theta[, 1])) * theta[, 1])
        P <- NQ \%*\% t(p)
        R <- as.matrix(rp)[, j] \%*\% t(JQ)
        gA <- pL * (R - P)/(P * (1 - P)) * dnorm(Z) * (TT)
        gB <- pL * (P - R)/(P * (1 - P)) * dnorm(Z)
        gA <- sum(rowSums(gA))
        gB <- sum(rowSums(gB))
        if (is.nan(gA) | is.nan(gB)) {
            print("Something was recorded as an Infinity...")
            print("p")
            print(p)
            print("A")
            print(aa)
            print("B")
            print(bb)
            print("Posterior")
            print(pL[1:10, ])
        }
        dL <- c(gA, gB)
    }
    if (tolower(settings$icc) == "logistic") {
        Z <- NQ \%*\% t(theta[, 1] - bb[j])
        if (settings$guess) {
            p <- cc[j] + (1 - cc[j]) * exp(1.7 * aa[j] * (theta[, 
                1] - bb[j]))/(1 + exp(1.7 * aa[j] * (theta[, 
                1] - bb[j])))
        }
        else {
            p <- exp(1.7 * aa[j] * (theta[, 1] - bb[j]))/(1 + 
                exp(1.7 * aa[j] * (theta[, 1] - bb[j])))
        }
        P <- NQ \%*\% t(p)
        R <- as.matrix(rp)[, j] \%*\% t(JQ)
        gA <- pL * (R - P) * (1.7) * Z
        gB <- pL * (R - P) * (-1.7) * aa[j]
        gA <- sum(rowSums(gA))
        gB <- sum(rowSums(gB))
        dL <- c(gA, gB)
        if (settings$guess) {
            gC <- pL * (1 - P)/(1 - cc[j])
            gC <- sum(rowSums(gC))
            dL <- c(dL, gC)
        }
    }
    return(dL)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
