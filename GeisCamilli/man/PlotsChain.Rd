\name{PlotsChain}
\alias{PlotsChain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PlotsChain(Aiter, Biter, settings, A, init)
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
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{init}{
%%     ~~Describe \code{init} here~~
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
function (Aiter, Biter, settings, A, init) 
{
    ifelse(grepl("\\.[Rr][Dd][Aa]", settings$simfile), load(file = settings$simfile), 
        load(file = paste(settings$simfile, ".rda", sep = "")))
    J <- nrow(gen.xi)
    N <- nrow(gen.rp)
    It <- ncol(Biter)
    par(mfrow = c(1, settings$Adim + 1), mar = c(3, 3, 1, 1))
    if (tolower(settings$fm) == "new") {
        out = princomp(covmat = gen.xi[, 1:settings$Adim] \%*\% 
            t(gen.xi[, 1:settings$Adim]), scores = FALSE, cor = FALSE)
        Avec = out$sd^2
        if (settings$Adim > 1) {
            Aload <- as.matrix(apply(out$loadings[, 1:settings$Adim], 
                2, function(x) (if (mean(x) > 0) {
                  return(x)
                }
                else {
                  return(-1 * x)
                })), nrow(gen.xi), settings$Adim)
            if (sum(Aload[, 2] * A[, 2]) < sum(Aload[, 2] * (-1) * 
                A[, 2])) 
                Aload[, 2] <- (-1) * Aload[, 2]
            Aplot = Aload[, 1:settings$Adim] \%*\% sqrt(diag(Avec[1:settings$Adim]))
        }
        else {
            Aload <- out$loadings[, 1]
            Aplot = Aload * sqrt(Avec[1])
        }
    }
    else if (tolower(settings$fm) \%in\% c("camilli", "eigen", 
        "pca")) {
        Aplot = A
    }
    for (i in 1:settings$Adim) {
        plot(Aiter[1, i, ], type = "n", xlim = c(0, length(Aiter[1, 
            1, ])), ylim = range(Aiter), main = paste("A", i))
        for (j in 1:J) {
            lines(Aiter[j, i, ], col = j)
        }
        if (tolower(settings$fm) == "new" | tolower(settings$fm) == 
            "camilli") {
            abline(h = gen.xi[, i])
        }
        else {
            if (settings$Adim == 1) {
                abline(h = Aplot)
            }
            else {
                abline(h = Aplot[, i])
            }
        }
        points(rep(1, J), init$XI[, i], pch = 16, cex = 1.5)
        abline(v = settings$burnin)
    }
    plot(Biter[1, ], type = "n", xlim = c(0, length(Biter[1, 
        ])), ylim = range(Biter), main = "B")
    for (j in 1:J) {
        lines(Biter[j, ], col = j)
    }
    if (!is.na(settings$simfile)) {
        abline(h = gen.xi[, ncol(gen.xi)])
    }
    points(rep(1, J), init$XI[, settings$Adim + 1], col = 1:J, 
        pch = 16, cex = 1.5)
    abline(v = settings$burnin)
    par(mfrow = c(2, settings$Adim + 2))
    SS <- 40
    window <- 30
    for (i in settings$Adim:1) {
        ASum <- apply(as.matrix(Aiter[, i, ]), 2, sum)[SS:settings$burnin]
        plot(SS:settings$burnin, ASum, type = "n", xlim = c(SS, 
            settings$burnin), ylim = range(ASum), main = paste("sum(A", 
            i, ")", sep = ""))
        lines(SS:settings$burnin, ASum)
    }
    BSum <- apply(as.matrix(abs(Biter)), 2, sum)[SS:settings$burnin]
    plot(SS:settings$burnin, BSum, type = "n", xlim = c(SS, settings$burnin), 
        ylim = range(BSum), main = "sum(abs(B))")
    lines(SS:settings$burnin, BSum)
    bSum <- apply(as.matrix(Biter), 2, sum)[SS:settings$burnin]
    plot(SS:settings$burnin, bSum, type = "n", xlim = c(SS, settings$burnin), 
        ylim = range(bSum), main = "sum(B)")
    lines(SS:settings$burnin, bSum)
    for (i in settings$Adim:1) {
        ASum <- apply(as.matrix(Aiter[, i, ]), 2, prod)[SS:settings$burnin]
        plot(SS:settings$burnin, ASum, type = "n", xlim = c(SS, 
            settings$burnin), ylim = range(ASum), main = paste("prod(A", 
            i, ")", sep = ""))
        lines(SS:settings$burnin, ASum)
    }
    BSum <- apply(as.matrix(abs(Biter)), 2, prod)[SS:settings$burnin]
    plot(SS:settings$burnin, BSum, type = "n", xlim = c(SS, settings$burnin), 
        ylim = range(BSum), main = "prod(abs(B))")
    lines(SS:settings$burnin, BSum)
    bSum <- apply(as.matrix(Biter), 2, prod)[SS:settings$burnin]
    plot(SS:settings$burnin, bSum, type = "n", xlim = c(SS, settings$burnin), 
        ylim = range(bSum), main = "prod(B)")
    lines(SS:settings$burnin, bSum)
    if (1/(It - settings$burnin)^settings$estgain < settings$eps) {
        VRB <- vector()
        VRb <- vector()
        VRA <- vector()
        VSB <- vector()
        VSb <- vector()
        VSA <- vector()
        for (j in 0:(settings$burnin - SS - 2 * window)) {
            VRB <- c(VRB, var(BSum[1 + SS + j - ceiling(window/2):0])/var(BSum[1 + 
                SS + j + ceiling(window/2):0]))
            VSB <- c(VSB, abs(sum(-1 * BSum[1 + SS + j - ceiling(window/2):0] + 
                BSum[1 + SS + j + ceiling(window/2):0])))
            VRb <- c(VRb, var(bSum[1 + SS + j - ceiling(window/2):0])/var(bSum[1 + 
                SS + j + ceiling(window/2):0]))
            VSb <- c(VSb, abs(sum(-1 * bSum[1 + SS + j - ceiling(window/2):0] + 
                bSum[1 + SS + j + ceiling(window/2):0])))
            VRA <- c(VRA, var(ASum[1 + SS + j - ceiling(window/2):0])/var(ASum[1 + 
                SS + j + ceiling(window/2):0]))
            VSA <- c(VSA, abs(sum(-1 * ASum[1 + SS + j - ceiling(window/2):0] + 
                ASum[1 + SS + j + ceiling(window/2):0])))
        }
        par(mfrow = c(2, 3))
        plot(SS + 0:(settings$burnin - SS - 2 * window), VRA, 
            main = "A1 Variance Ratio")
        plot(SS + 0:(settings$burnin - SS - 2 * window), VRB, 
            main = "B Variance Ratio")
        plot(SS + 0:(settings$burnin - SS - 2 * window), VRb, 
            main = "abs(B) Variance Ratio")
        plot(SS + 0:(settings$burnin - SS - 2 * window), VSA, 
            main = "A1 Signal Search")
        plot(SS + 0:(settings$burnin - SS - 2 * window), VSB, 
            main = "B Signal Search")
        plot(SS + 0:(settings$burnin - SS - 2 * window), VSb, 
            main = "abs(B) Signal Search")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
