\name{PlotRecon}
\alias{PlotRecon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PlotRecon(gen, fit, estfile, simfile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gen}{
%%     ~~Describe \code{gen} here~~
}
  \item{fit}{
%%     ~~Describe \code{fit} here~~
}
  \item{estfile}{
%%     ~~Describe \code{estfile} here~~
}
  \item{simfile}{
%%     ~~Describe \code{simfile} here~~
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
function (gen, fit, estfile, simfile) 
{
    Err <- abs(fit$xiError)^0.5
    print(Err)
    BErr <- Err[, ncol(gen$XI)]
    AErr <- Err[, 1:(ncol(gen$XI) - 1)]
    if (!is.na(estfile) & !is.na(simfile)) {
        if (grepl("\\.[Rr][Dd][Aa]", estfile)) {
            efilename = estfile
        }
        else {
            efilename = paste(estfile, ".rda", sep = "")
        }
        load(file = efilename)
        if (grepl("\\.[Rr][Dd][Aa]", simfile)) {
            sfilename = simfile
        }
        else {
            sfilename = paste(simfile, ".rda", sep = "")
        }
        load(file = sfilename)
        par(mfrow = c(1, ncol(gen$XI)), mar = c(3, 3, 1, 1))
        for (q in 1:(ncol(gen$XI) - 1)) {
            plot(MCMCDATA$Aiter[1, q, ], type = "n", main = "A Parameters", 
                xlab = "Iterations", ylab = "Value", ylim = c(-0.3, 
                  0.3) + range(gen$XI[, q]))
            for (i in 1:nrow(gen$XI)) {
                lines(MCMCDATA$Aiter[i, q, ], col = i)
                abline(h = gen$XI[i, q], col = i, lwd = 1.5)
            }
        }
        plot(MCMCDATA$Biter[1, ], type = "n", main = "B Parameters", 
            xlab = "Iterations", ylab = "Value", ylim = 1.1 * 
                range(gen$XI[, ncol(gen$XI)]))
        for (i in 1:nrow(gen$XI)) {
            lines(MCMCDATA$Biter[i, ], col = i)
            abline(h = gen$XI[i, ncol(gen$XI)], col = i, lwd = 1.5)
        }
        if (ncol(gen$XI) > 2) {
            par(mfrow = c(2, 2))
            plot(MCMCDATA$LLiter, type = "n", main = "Likelihood over iterations")
            lines(MCMCDATA$LLiter)
            VA <- colSums((as.matrix(MCMCDATA$Aiter[, 1, ]) - 
                as.matrix(gen$XI[, 1]) \%*\% rep(1, length(MCMCDATA$Aiter[1, 
                  1, ])))^2)
            for (q in 2:(ncol(gen$XI) - 1)) VA <- VA + colSums((as.matrix(MCMCDATA$Aiter[, 
                q, ]) - as.matrix(gen$XI[, q]) \%*\% rep(1, length(MCMCDATA$Aiter[1, 
                q, ])))^2)
            VB <- colSums((MCMCDATA$Biter - as.matrix(gen$XI[, 
                2]) \%*\% rep(1, ncol(MCMCDATA$Biter)))^2)
            plot(VA, VB, main = "Total SSE, B vs. A by Iteration", 
                type = "n")
            lines(VA, VB)
            points(VA[length(VA)], VB[length(VB)], pch = 19, 
                col = "red")
            plot(VA[-1], MCMCDATA$LLiter, type = "n", main = "Likelihood vs SSE of Parameter A")
            lines(VA[-1], MCMCDATA$LLiter)
            points(VA[length(VA)], MCMCDATA$LLiter[length(MCMCDATA$LLiter)], 
                col = "red", pch = 19)
            plot(VB[-1], MCMCDATA$LLiter, type = "n", main = "Likelihood vs SSE of Parameter B")
            lines(VB[-1], MCMCDATA$LLiter)
            points(VB[length(VB)], MCMCDATA$LLiter[length(MCMCDATA$LLiter)], 
                col = "red", pch = 19)
        }
    }
    par(mfrow = c(ncol(gen$XI), 1), mar = c(4, 4, 3, 2))
    XL <- ""
    for (i in 1:(ncol(gen$XI) - 1)) XL <- paste(XL, "A", i, "(1...J);", 
        sep = "")
    plot(as.vector(fit$xi) - as.vector(gen$XI), main = "Parameter Estimate Differences", 
        ylab = "XI_hat - XI_gen", xlab = paste(XL, "B(1...J)", 
            sep = ""), ylim = range(c(as.vector(fit$xi) - as.vector(gen$XI) + 
            2 * c(as.vector(AErr), as.vector(BErr)), as.vector(fit$xi) - 
            as.vector(gen$XI) - 2 * c(as.vector(AErr), as.vector(BErr)))), 
        col = "red")
    arrows(1:length(as.vector(fit$xi)), as.vector(fit$xi) - as.vector(gen$XI) - 
        2 * c(as.vector(AErr), as.vector(BErr)), 1:length(as.vector(fit$xi)), 
        as.vector(fit$xi) - as.vector(gen$XI) + 2 * c(as.vector(AErr), 
            as.vector(BErr)), code = 3, angle = 90, length = 0.07)
    abline(h = 0)
    if (max(abs(as.vector(fit$xi) - as.vector(gen$XI))) > 0.3) {
        text(which(abs(as.vector(fit$xi) - as.vector(gen$XI)) > 
            0.3), (as.vector(fit$xi) - as.vector(gen$XI))[which(abs(as.vector(fit$xi) - 
            as.vector(gen$XI)) > 0.3)] + 0.05, ((which(abs(as.vector(fit$xi) - 
            as.vector(gen$XI)) > 0.3) - 1)\%\%nrow(gen$XI)) + 1)
    }
    axes <- c(paste("A", 1:(ncol(gen$XI) - 1)), "B")
    for (i in 1:(ncol(gen$XI) - 1)) {
        plot(gen$XI[, i], gen$XI[, i + 1], main = "Target is the Black Circle, Fit is Red, Line is Bias", 
            xlab = axes[i], ylab = axes[i + 1], pch = 19, ylim = 1.15 * 
                range(fit$xi[, i + 1]), xlim = c(0.5, 1.15) * 
                range(fit$xi[, i]))
        arrows(gen$XI[, i], gen$XI[, i + 1], fit$xi[, i], gen$XI[, 
            i + 1], length = 0)
        arrows(gen$XI[, i], gen$XI[, i + 1], gen$XI[, i], fit$xi[, 
            i + 1], length = 0)
        arrows(fit$xi[, i], fit$xi[, i + 1], fit$xi[, i], gen$XI[, 
            i + 1], length = 0)
        arrows(fit$xi[, i], fit$xi[, i + 1], gen$XI[, i], fit$xi[, 
            i + 1], length = 0)
        points(fit$xi[, i], fit$xi[, i + 1], pch = 19, col = "red")
        points(fit$xi[, i], fit$xi[, i + 1], cex = 1.1, col = "red")
    }
    if (ncol(gen$XI) > 2) {
        par(mfrow = c(ncol(gen$XI) - 1, 10))
        for (i in 1:(ncol(gen$XI) - 1)) {
            plot(gen$THETA[, i], fit$TROT[, i], main = paste("Theta", 
                i, ": Boot"), xlab = "Theta_gen", ylab = paste("Theta Fit (Bootstrapped Z,Theta)", 
                sep = ""))
            abline(a = 0, b = 1)
            TSD <- format(2 * sd(fit$TROT[, i] - gen$THETA[, 
                i]), digits = 4)
            hist(fit$TROT[, i] - gen$THETA[, i], breaks = 20, 
                main = paste("Boot Theta", i, "- Theta_gen\n2*SE =", 
                  TSD), xlab = "Theta Fit - Theta_gen")
            plot(gen$THETA[, i], fit$TMAPROT[, i], main = paste("Theta", 
                i, ": MAP(pca)"), xlab = "Theta_gen", ylab = paste("Theta Fit (MAP(pca))", 
                sep = ""))
            abline(a = 0, b = 1)
            TSD <- format(2 * sd(fit$TMAPROT[, i] - gen$THETA[, 
                i]), digits = 4)
            hist(fit$TMAPROT[, i] - gen$THETA[, i], breaks = 20, 
                main = paste("MAP(pca) Theta", i, "- Theta_gen\n2*SE =", 
                  TSD), xlab = "Theta Fit - Theta_gen")
            plot(gen$THETA[, i], fit$TRMAP[, i], main = paste("Theta", 
                i, ": MAP(Rot)"), xlab = "Theta_gen", ylab = paste("Theta Fit (MAP(Rot))", 
                sep = ""))
            abline(a = 0, b = 1)
            TSD <- format(2 * sd(fit$TRMAP[, i] - gen$THETA[, 
                i]), digits = 4)
            hist(fit$TRMAP[, i] - gen$THETA[, i], breaks = 20, 
                main = paste("MAP(Rot) Theta", i, "- Theta_gen\n2*SE =", 
                  TSD), xlab = "Theta Fit - Theta_gen")
            plot(gen$THETA[, i], fit$ThetaFix$FT[, i], main = paste("Theta", 
                i, ": Fixed"), xlab = "Theta_gen", ylab = "Theta Fit, Fixed Parameters")
            abline(a = 0, b = 1)
            TSD <- format(2 * sd(fit$ThetaFix$FT[, i] - gen$THETA[, 
                i]), digits = 4)
            MSD <- format(2 * mean(fit$ThetaFix$FTErr[, i]), 
                digits = 4)
            hist(fit$ThetaFix$FT[, i] - gen$THETA[, i], breaks = 20, 
                main = paste("Fixed P Theta", i, "- Theta_gen\n2*SE =", 
                  TSD, "; Iterative 2*SE =", MSD), xlab = "Theta Fit - Theta_gen")
            plot(gen$THETA[, i], fit$EmpSE$MET[, i], main = paste("Theta", 
                i, ": Empirical"), xlab = "Theta_gen", ylab = "Theta Fit, Empirical")
            abline(a = 0, b = 1)
            TSD <- format(2 * sd(fit$EmpSE$MET[, i] - gen$THETA[, 
                i]), digits = 4)
            MSD <- format(2 * mean(fit$EmpSE$SET[, i]), digits = 4)
            hist(fit$EmpSE$MET[, i] - gen$THETA[, i], breaks = 20, 
                main = paste("Empirical Theta", i, "- Theta_gen\n2*SE =", 
                  TSD, "; Iterative 2*SE =", MSD), xlab = "Theta Fit - Theta_gen")
        }
    }
    else {
        par(mfrow = c(ncol(gen$XI) - 1, 6))
        plot(gen$THETA, fit$THAT[, 1], main = "Theta (Boot) vs. Theta_gen", 
            xlab = "Theta_gen", ylab = paste("Theta Fit (Boot)", 
                sep = ""))
        abline(a = 0, b = 1)
        TSD <- format(2 * sd(fit$THAT[, 1] - gen$THETA), digits = 4)
        hist(fit$THAT[, 1] - gen$THETA, breaks = 20, main = paste("Theta-Theta_gen\n2*SE =", 
            TSD), xlab = "Theta Fit - Theta_gen")
        plot(gen$THETA, fit$ThetaFix$FT, main = "Fixed Param Theta vs. Theta_gen", 
            xlab = "Theta_gen", ylab = "Theta Fit (Fixed P)")
        abline(a = 0, b = 1)
        TSD <- format(2 * sd(fit$ThetaFix$FT - gen$THETA), digits = 4)
        MSD <- format(2 * mean(fit$ThetaFix$FTErr), digits = 4)
        hist(fit$ThetaFix$FT - gen$THETA, breaks = 20, main = paste("Theta-Theta_gen\n2*SE =", 
            TSD, "; Iterative 2*SE =", MSD), xlab = "Theta Fit - Theta_gen")
        plot(gen$THETA, fit$EmpSE$MET, main = "Floating Param Theta vs. Theta_gen", 
            xlab = "Theta_gen", ylab = "Theta Fit (Float Params)")
        abline(a = 0, b = 1)
        TSD <- format(2 * sd(fit$EmpSE$MET - gen$THETA), digits = 4)
        MSD <- format(2 * mean(fit$EmpSE$SET), digits = 4)
        hist(fit$EmpSE$MET - gen$THETA, breaks = 20, main = paste("Theta-Theta_gen\nSE =", 
            TSD, "; Iterative 2*SE =", MSD), xlab = "Theta Fit - Theta_gen")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
