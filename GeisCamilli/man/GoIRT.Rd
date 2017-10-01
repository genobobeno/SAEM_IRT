\name{GoIRT}
\alias{GoIRT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GoIRT(rp, init = Init, settings = settings)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rp}{
%%     ~~Describe \code{rp} here~~
}
  \item{init}{
%%     ~~Describe \code{init} here~~
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
function (rp, init = Init, settings = settings) 
{
    stopifnot(length(settings$tmu) == settings$Adim, ncol(init$XI) == 
        (1 + settings$Adim + (settings$guess + 0)))
    prior <- list(tmu = settings$tmu, tsigma = settings$tsigma)
    J <- ncol(rp)
    N <- nrow(rp)
    THat <- init$THat
    xi0 <- init$XI
    Theta = GetQuad(settings)
    print(paste("Students in Sample:", N))
    Q <- length(Theta)
    it = 1
    repeat {
        print(it)
        FL <- IRTFullLL(xi = xi0, rp = rp, settings = settings)
        PT <- Posterior(L = FL, xi = xi0, rp = rp, settings = settings)
        errors <- mat.or.vec(J, ncol(xi0))
        xiN <- xiI <- xi0
        for (j in 1:J) {
            print(paste("Item", j))
            ncycle <- 1
            repeat {
                grad <- Gradient(j, xi = xiN, rp = rp, pL = PT, 
                  settings = settings)
                hess <- Hessian(j, xi = xiN, rp = rp, pL = PT, 
                  settings = settings)
                if (settings$guess) {
                  info <- matrix(c(hess[1], hess[4], hess[5], 
                    hess[4], hess[2], hess[6], hess[5], hess[6], 
                    hess[3]), nrow = 3, ncol = 3)
                  if (is.nan(info[3, 3])) 
                    info[3, 3] <- 1
                  check <- ginv(info) \%*\% grad
                }
                else {
                  info <- matrix(c(hess[1], hess[3], hess[3], 
                    hess[2]), nrow = length(grad), ncol = length(grad))
                  print(info)
                  check <- ginv(info) \%*\% grad
                }
                xiN[j, ] <- xiI[j, ] - settings$estgain * check
                print(paste("Newton difference:", abs(max(xiI[j, 
                  ] - xiN[j, ])), "after", ncycle, "Newton Cycles"))
                if (abs(max(xiN[j, ] - xiI[j, ])) < settings$eps) {
                  Err <- (-1) * ginv(info)
                  if (settings$guess) {
                    errors[j, ] <- c(Err[1, 1], Err[2, 2], Err[3, 
                      3])
                  }
                  else {
                    errors[j, ] <- c(Err[1, 1], Err[2, 2])
                  }
                  if (It\%\%10 == 0 && settings$plots) 
                    plot(exp(1.7 * xiN[j, 1] * (Theta[, 1] - 
                      xiN[j, 2]))/(1 + exp(1.7 * xiN[j, 1] * 
                      (Theta[, 1] - xiN[j, 2]))) ~ Theta[, 1], 
                      main = paste("Item", j), ylab = "Probability")
                  xiI[j, 1] <- xiN[j, 1]
                  xiI[j, 2] <- xiN[j, 2]
                  if (settings$guess) 
                    xiI[j, 3] <- xiN[j, 3]
                  break
                }
                if (xiN[j, 1] < 0.1) {
                  xiN[j, 1] <- runif(1, min = 0.25, max = 0.7)
                }
                if (xiN[j, 1] > 2.8) {
                  xiN[j, 1] <- runif(1, min = 1, max = 1.7)
                }
                if (xiN[j, 2] < -3.5) {
                  xiN[j, 2] <- runif(1, min = -3, max = -1)
                }
                if (xiN[j, 2] > 3.5) {
                  xiN[j, 2] <- runif(1, min = 1, max = 3)
                }
                if (settings$guess) 
                  if (xiN[j, 3] > 0.4 | xiN[j, 3] < 0) 
                    xiN[j, 3] <- runif(1, min = 0.05, max = 0.15)
                xiI <- xiN
                ncycle <- ncycle + 1
            }
        }
        print(paste("Max dA:", max(abs(xi0[, 1] - xiI[, 1])), 
            "   Max dB:", max(abs(xi0[, 2] - xiI[, 2]))))
        if (max(sqrt((xi0 - xiI)^2)) < settings$eps) {
            xi0 <- xiI
            EAP <- PT \%*\% Theta[, 1]
            THat <- EAP[, 1]
            print(paste("After", it, "iterations, we have converged"))
            break
        }
        xi0 <- xiI
        it <- it + 1
    }
    FitDATA <- list(xi = xi0, xiError = errors, THAT = THat)
    if (grepl("\\.[Rr][Dd][Aa]", settings$estfile)) {
        filename = settings$estfile
    }
    else {
        filename = paste(settings$estfile, ".rda", sep = "")
    }
    save(FitDATA, settings, file = filename)
    return(FitDATA)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
