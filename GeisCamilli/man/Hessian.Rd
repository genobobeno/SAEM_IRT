\name{Hessian}
\alias{Hessian}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Hessian(j, xi, rp, pL, settings)
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
        dA <- dnorm(Z) * (TT)
        dB <- (-1) * dnorm(Z)
        dAdB <- Z * TT * dnorm(Z)
        ddA <- (-1) * Z^2 * TT^2 * dnorm(Z)
        ddB <- Z * dnorm(Z)
        ddLA <- pL/(P * (1 - P)) * ((R - P) * ddA - ((1 - 2 * 
            P) * (R - P)/(P * (1 - P)) + 1) * dA^2)
        ddLB <- pL/(P * (1 - P)) * ((R - P) * ddB - ((1 - 2 * 
            P) * (R - P)/(P * (1 - P)) + 1) * dB^2)
        ddLAB <- pL/(P * (1 - P)) * ((R - P) * dAdB - ((1 - 2 * 
            P) * (R - P)/(P * (1 - P)) + 1) * dB * dA)
        ddLA <- sum(rowSums(ddLA))
        ddLB <- sum(rowSums(ddLB))
        ddLAB <- sum(rowSums(ddLAB))
        ddL <- c(ddLA, ddLB, ddLAB)
    }
    else if (tolower(settings$icc) == "logistic") {
        Z <- NQ \%*\% t(theta[, 1] - bb[j])
        p <- exp(1.7 * aa[j] * (theta[, 1] - bb[j]))/(1 + exp(1.7 * 
            aa[j] * (theta[, 1] - bb[j])))
        P <- NQ \%*\% t(p)
        R <- as.matrix(rp)[, j] \%*\% t(JQ)
        if (settings$guess) {
            gA <- (1.7/(1 - cc[j]))^2 * pL * (R * cc[j]/P^2 - 
                1) * Z^2 * (P - cc[j]) * (1 - P)
            gB <- (1.7 * aa[j]/(1 - cc[j]))^2 * pL * (R * cc[j]/P^2 - 
                1) * (P - cc[j]) * (1 - P)
            gC <- pL/(1 - cc[j])^2 * (R/P - 1) * (R * (1 - P)/P^2)
            gAB <- (-1.7)/(1 - cc[j]) * pL * (P - cc[j]) * ((R/P - 
                1) + 1.7 * aa[j]/(1 - cc[j]) * Z * (1 - P)/P * 
                (R * cc[j]/P - P))
            gAC <- (-1.7)/(1 - cc[j])^2 * pL * Z * (P - cc[j]) * 
                (1 - P)/P^2 * R
            gBC <- 1.7 * aa[j]/(1 - cc[j])^2 * pL * (P - cc[j]) * 
                (1 - P)/P^2 * R
            ddA <- sum(rowSums(gA))
            ddB <- sum(rowSums(gB))
            ddC <- sum(rowSums(gC))
            dAdB <- sum(rowSums(gAB))
            dAdC <- sum(rowSums(gAC))
            dBdC <- sum(rowSums(gBC))
            ddL <- c(ddA, ddB, ddC, dAdB, dAdC, dBdC)
        }
        else {
            gA <- pL * Z^2 * P * (1 - P)
            ddA <- (-1.7 * 1.7) * sum(rowSums(gA))
            gB <- pL * P * (1 - P) * aa[j]^2
            ddB <- (-1.7 * 1.7) * sum(rowSums(gB))
            gAB <- pL * ((R - P) - 1.7 * aa[j] * Z * P * (1 - 
                P))
            dAdB <- (-1.7) * sum(rowSums(gAB))
            ddL <- c(ddA, ddB, dAdB)
        }
    }
    return(ddL)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
