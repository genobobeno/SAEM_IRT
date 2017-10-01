\name{PlotsCCF}
\alias{PlotsCCF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PlotsCCF(Aiter, Biter, settings, ItCC, ItAC)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Aiter}{
%%     ~~Describe \code{Aiter} here~~
}
  \item{Biter}{
%%     ~~Describe \code{Biter} here~~
}
  \item{settings}{
%%     ~~Describe \code{settings} here~~
}
  \item{ItCC}{
%%     ~~Describe \code{ItCC} here~~
}
  \item{ItAC}{
%%     ~~Describe \code{ItAC} here~~
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
function (Aiter, Biter, settings, ItCC, ItAC) 
{
    J <- nrow(Biter)
    It <- ncol(Biter)
    start <- 1
    par(mfrow = c(floor(J/5), floor(J/5)), mar = c(2, 2, 0.5, 
        0.5))
    TotCC = 0
    TotAC = rep(1, 10)
    for (j in 1:J) {
        for (jj in 1:J) {
            if (j != jj) {
                A1cf <- ccf(Aiter[j, 1, start:It], Aiter[jj, 
                  1, start:It], lag.max = 10, cex = 0.8, ylim = c(-1, 
                  1))
                text(-5, -0.5, paste("A1:", j, "-A1:", jj, sep = ""))
                TotCC = TotCC + max(abs(A1cf$acf))
            }
            else {
                A1cf <- ccf(Aiter[j, 1, start:It], Aiter[jj, 
                  2, start:It], lag.max = 10, cex = 0.8, ylim = c(-1, 
                  1))
                text(-5, -0.5, paste("A1:", j, "-A2:", jj, sep = ""))
                TotCC = TotCC + max(abs(A1cf$acf))
                A1cf <- acf(Aiter[jj, 1, 1:It], plot = FALSE, 
                  cex = 0.8, ylim = c(-1, 1), type = "correlation")
                TotAC = TotAC * abs(A1cf$acf[2:11])
            }
        }
    }
    for (j in 1:J) {
        for (jj in 1:J) {
            if (j != jj) {
                A2cf <- ccf(Aiter[j, 2, start:It], Aiter[jj, 
                  2, start:It], lag.max = 10, cex = 0.8, ylim = c(-1, 
                  1))
                text(-5, -0.5, paste("A2:", j, "-A2:", jj, sep = ""))
                TotCC = TotCC + max(abs(A2cf$acf))
            }
            else {
                A2cf <- ccf(Aiter[j, 2, start:It], Aiter[jj, 
                  1, start:It], lag.max = 10, cex = 0.8, ylim = c(-1, 
                  1))
                text(-5, -0.5, paste("A2:", j, "-A1:", jj, sep = ""))
                TotCC = TotCC + max(abs(A2cf$acf))
                A2cf <- acf(Aiter[jj, 2, 1:It], plot = FALSE, 
                  cex = 0.8, ylim = c(-1, 1), type = "correlation")
                TotAC = TotAC * abs(A2cf$acf[2:11])
            }
        }
    }
    for (j in 1:J) {
        for (jj in 1:J) {
            A1cf <- ccf(Aiter[j, 1, start:It], Biter[jj, start:It], 
                lag.max = 10, cex = 0.8, ylim = c(-1, 1))
            text(-5, -0.5, paste("A1:", j, "-B:", jj, sep = ""))
            TotCC = TotCC + max(abs(A1cf$acf))
        }
    }
    for (j in 1:J) {
        for (jj in 1:J) {
            A2cf <- ccf(Aiter[j, 2, start:It], Biter[jj, start:It], 
                lag.max = 10, cex = 0.8, ylim = c(-1, 1))
            text(-5, -0.5, paste("A2:", j, "-B:", jj, sep = ""))
            TotCC = TotCC + max(abs(A2cf$acf))
        }
    }
    for (j in 1:J) {
        for (jj in 1:J) {
            if (j != jj) {
                Bcf <- ccf(Biter[j, start:It], Biter[jj, start:It], 
                  lag.max = 10, cex = 0.8, ylim = c(-1, 1))
                text(-5, -0.5, paste("B:", j, "-B:", jj, sep = ""))
                TotCC = TotCC + max(abs(Bcf$acf))
            }
            else {
                Bcf <- acf(Biter[jj, start:It], cex = 0.8, ylim = c(-1, 
                  1), type = "correlation")
                text(-5, -0.5, paste("ACF B:", jj, sep = ""))
                TotAC = TotAC * abs(Bcf$acf[2:11])
            }
        }
    }
    ItCC <<- c(ItCC, TotCC)
    ItAC <<- c(ItAC, sum(TotAC))
    if (It > 50) {
        par(mfrow = c(2, 1))
        plot(ItCC, main = "Total Cross Correlation starting from Iteration 30")
        plot(ItAC, main = "Total Auto Correlation starting from Iteration 30")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
