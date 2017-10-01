\name{PlotErrors}
\alias{PlotErrors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PlotErrors(Fit, setting = settings)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Fit}{
%%     ~~Describe \code{Fit} here~~
}
  \item{setting}{
%%     ~~Describe \code{setting} here~~
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
function (Fit, setting = settings) 
{
    load("RMSE.rda")
    items = length(Fit$B)
    examinees = nrow(as.matrix(Fit$THAT))
    if (setting$Adim == 1) {
        par(mfrow = c(2, 2))
        plot(Fit$A, Fit$EmpSE$MCSA, xlab = "A (Generated)", ylab = "Standard Error", 
            pch = 19, main = "dA: Comparing initseq() SE to RMSE")
        abline(lm(Fit$EmpSE$MCSA ~ Fit$A)$coefficients, lwd = 1.5)
        points(RMSE$stats[RMSE$stats$N == examinees & RMSE$stats$J == 
            items, c("A1", "A1_sd")], col = 2, pch = 19)
        abline(lm(RMSE$stats[RMSE$stats$N == examinees & RMSE$stats$J == 
            items, c("A1_sd", "A1")])$coefficients, lwd = 1.5, 
            col = 2)
        legend("topleft", c("MCMC", "RMSE"), pch = 19, col = c(1, 
            2))
        plot(Fit$B, Fit$EmpSE$MCSB, xlab = "B (Generated)", ylab = "Standard Error", 
            pch = 19, main = "dB: Comparing initseq() SE to RMSE")
        lines(Fit$B[order(Fit$B)], predict(lm(Fit$EmpSE$MCSB ~ 
            Fit$B + I(Fit$B^2)))[order(Fit$B)], lwd = 1.5)
        points(RMSE$stats[RMSE$stats$N == examinees & RMSE$stats$J == 
            items, c("B", "B_sd")], col = 2, pch = 19)
        lines(RMSE$stats[RMSE$stats$N == examinees & RMSE$stats$J == 
            items, "B"][order(RMSE$stats[RMSE$stats$N == examinees & 
            RMSE$stats$J == items, "B"])], predict(lm(RMSE$stats[RMSE$stats$N == 
            examinees & RMSE$stats$J == items, "B_sd"] ~ RMSE$stats[RMSE$stats$N == 
            examinees & RMSE$stats$J == items, "B"] + I(RMSE$stats[RMSE$stats$N == 
            examinees & RMSE$stats$J == items, "B"]^2)))[order(RMSE$stats[RMSE$stats$N == 
            examinees & RMSE$stats$J == items, "B"])], lwd = 1.5, 
            col = 2)
        legend("topleft", c("MCMC", "RMSE"), pch = 19, col = c(1, 
            2))
        plot(RMSE$stats[RMSE$stats$N == examinees & RMSE$stats$J == 
            items, "A1_sd"], Fit$EmpSE$MCSA, xlab = "RMSE (A)", 
            ylab = "MCMC Standard Error", pch = 19, main = "A: MCMC SE vs. RMSE")
        abline(a = 0, b = 1)
        plot(RMSE$stats[RMSE$stats$N == examinees & RMSE$stats$J == 
            items, "B_sd"], Fit$EmpSE$MCSB, xlab = "RMSE (B)", 
            ylab = "MCMC Standard Error", pch = 19, main = "B: MCMC SE vs. RMSE")
        abline(a = 0, b = 1)
        par(mfrow = c(2, 2))
        plot(Fit$A, sqrt(Fit$EmpSE$VARLMIo[, 1]), xlab = "A (Generated)", 
            ylab = "Standard Error", pch = 19, main = "dA: Comparing Ogive Emp SE to RMSE")
        abline(lm(sqrt(Fit$EmpSE$VARLMIo[, 1]) ~ Fit$A)$coefficients, 
            lwd = 1.5)
        points(RMSE$stats[RMSE$stats$N == examinees & RMSE$stats$J == 
            items, c("A1", "A1_sd")], col = 2, pch = 19)
        abline(lm(RMSE$stats[RMSE$stats$N == examinees & RMSE$stats$J == 
            items, c("A1_sd", "A1")])$coefficients, lwd = 1.5, 
            col = 2)
        legend("topleft", c("Ogive Emp SE", "RMSE"), pch = 19, 
            col = c(1, 2))
        plot(Fit$B, sqrt(abs(Fit$EmpSE$VARLMIo[, 2])), xlab = "B (Generated)", 
            ylab = "Standard Error", pch = 19, main = "dB: Comparing Ogive Emp SE to RMSE")
        lines(Fit$B[order(Fit$B)], predict(lm(sqrt(abs(Fit$EmpSE$VARLMIo[, 
            2])) ~ Fit$B + I(Fit$B^2)))[order(Fit$B)], lwd = 1.5)
        points(RMSE$stats[RMSE$stats$N == examinees & RMSE$stats$J == 
            items, c("B", "B_sd")], col = 2, pch = 19)
        lines(RMSE$stats[RMSE$stats$N == examinees & RMSE$stats$J == 
            items, "B"][order(RMSE$stats[RMSE$stats$N == examinees & 
            RMSE$stats$J == items, "B"])], predict(lm(RMSE$stats[RMSE$stats$N == 
            examinees & RMSE$stats$J == items, "B_sd"] ~ RMSE$stats[RMSE$stats$N == 
            examinees & RMSE$stats$J == items, "B"] + I(RMSE$stats[RMSE$stats$N == 
            examinees & RMSE$stats$J == items, "B"]^2)))[order(RMSE$stats[RMSE$stats$N == 
            examinees & RMSE$stats$J == items, "B"])], lwd = 1.5, 
            col = 2)
        legend("topleft", c("Ogive Emp SE", "RMSE"), pch = 19, 
            col = c(1, 2))
        plot(RMSE$stats[RMSE$stats$N == examinees & RMSE$stats$J == 
            items, "A1_sd"], sqrt(Fit$EmpSE$VARLMIo[, 1]), xlab = "RMSE (A)", 
            ylab = "Ogive Emp Standard Error", pch = 19, main = "A: Ogive Emp SE vs. RMSE")
        abline(a = 0, b = 1)
        plot(RMSE$stats[RMSE$stats$N == examinees & RMSE$stats$J == 
            items, "B_sd"], sqrt(Fit$EmpSE$VARLMIo[, 2]), xlab = "RMSE (B)", 
            ylab = "Ogive Emp Standard Error", pch = 19, main = "B: Ogive Emp SE vs. RMSE")
        abline(a = 0, b = 1)
        par(mfrow = c(2, 2))
        plot(Fit$A, sqrt(Fit$EmpSE$VARLMIi[, 1]), xlab = "A (Generated)", 
            ylab = "Standard Error", pch = 19, main = "dA: Comparing IRT Emp SE to RMSE")
        abline(lm(sqrt(Fit$EmpSE$VARLMIi[, 1]) ~ Fit$A)$coefficients, 
            lwd = 1.5)
        points(RMSE$stats[RMSE$stats$N == examinees & RMSE$stats$J == 
            items, c("A1", "A1_sd")], col = 2, pch = 19)
        abline(lm(RMSE$stats[RMSE$stats$N == examinees & RMSE$stats$J == 
            items, c("A1_sd", "A1")])$coefficients, lwd = 1.5, 
            col = 2)
        legend("topleft", c("IRT Emp SE", "RMSE"), pch = 19, 
            col = c(1, 2))
        plot(Fit$B, sqrt(abs(Fit$EmpSE$VARLMIi[, 2])), xlab = "B (Generated)", 
            ylab = "Standard Error", pch = 19, main = "dB: Comparing IRT Emp SE to RMSE")
        lines(Fit$B[order(Fit$B)], predict(lm(sqrt(abs(Fit$EmpSE$VARLMIi[, 
            2])) ~ Fit$B + I(Fit$B^2)))[order(Fit$B)], lwd = 1.5)
        points(RMSE$stats[RMSE$stats$N == examinees & RMSE$stats$J == 
            items, c("B", "B_sd")], col = 2, pch = 19)
        lines(RMSE$stats[RMSE$stats$N == examinees & RMSE$stats$J == 
            items, "B"][order(RMSE$stats[RMSE$stats$N == examinees & 
            RMSE$stats$J == items, "B"])], predict(lm(RMSE$stats[RMSE$stats$N == 
            examinees & RMSE$stats$J == items, "B_sd"] ~ RMSE$stats[RMSE$stats$N == 
            examinees & RMSE$stats$J == items, "B"] + I(RMSE$stats[RMSE$stats$N == 
            examinees & RMSE$stats$J == items, "B"]^2)))[order(RMSE$stats[RMSE$stats$N == 
            examinees & RMSE$stats$J == items, "B"])], lwd = 1.5, 
            col = 2)
        legend("topleft", c("IRT Emp SE", "RMSE"), pch = 19, 
            col = c(1, 2))
        plot(RMSE$stats[RMSE$stats$N == examinees & RMSE$stats$J == 
            items, "A1_sd"], sqrt(Fit$EmpSE$VARLMIi[, 1]), xlab = "RMSE (A)", 
            ylab = "IRT Emp Standard Error", pch = 19, main = "A: IRT Emp SE vs. RMSE")
        abline(a = 0, b = 1)
        plot(RMSE$stats[RMSE$stats$N == examinees & RMSE$stats$J == 
            items, "B_sd"], sqrt(Fit$EmpSE$VARLMIi[, 2]), xlab = "RMSE (B)", 
            ylab = "IRT Emp Standard Error", pch = 19, main = "B: IRT Emp SE vs. RMSE")
        abline(a = 0, b = 1)
        par(mfrow = c(2, 2))
        plot(Fit$A, sqrt(abs(Fit$xiError[, 1])), xlab = "A (Generated)", 
            ylab = "Standard Error", pch = 19, main = "dA: Comparing IRT SE to RMSE")
        abline(lm(sqrt(abs(Fit$xiError[, 1])) ~ Fit$A)$coefficients, 
            lwd = 1.5)
        points(RMSE$stats[RMSE$stats$N == examinees & RMSE$stats$J == 
            items, c("A1", "A1_sd")], col = 2, pch = 19)
        abline(lm(RMSE$stats[RMSE$stats$N == examinees & RMSE$stats$J == 
            items, c("A1_sd", "A1")])$coefficients, lwd = 1.5, 
            col = 2)
        legend("topleft", c("IRT SE", "RMSE"), pch = 19, col = c(1, 
            2))
        plot(Fit$B, sqrt(abs(Fit$xiError[, 2])), xlab = "B (Generated)", 
            ylab = "Standard Error", pch = 19, main = "dB: Comparing IRT SE to RMSE")
        lines(Fit$B[order(Fit$B)], predict(lm(sqrt(abs(Fit$xiError[, 
            2])) ~ Fit$B + I(Fit$B^2)))[order(Fit$B)], lwd = 1.5)
        points(RMSE$stats[RMSE$stats$N == examinees & RMSE$stats$J == 
            items, c("B", "B_sd")], col = 2, pch = 19)
        lines(RMSE$stats[RMSE$stats$N == examinees & RMSE$stats$J == 
            items, "B"][order(RMSE$stats[RMSE$stats$N == examinees & 
            RMSE$stats$J == items, "B"])], predict(lm(RMSE$stats[RMSE$stats$N == 
            examinees & RMSE$stats$J == items, "B_sd"] ~ RMSE$stats[RMSE$stats$N == 
            examinees & RMSE$stats$J == items, "B"] + I(RMSE$stats[RMSE$stats$N == 
            examinees & RMSE$stats$J == items, "B"]^2)))[order(RMSE$stats[RMSE$stats$N == 
            examinees & RMSE$stats$J == items, "B"])], lwd = 1.5, 
            col = 2)
        legend("topleft", c("IRT SE", "RMSE"), pch = 19, col = c(1, 
            2))
        plot(RMSE$stats[RMSE$stats$N == examinees & RMSE$stats$J == 
            items, "A1_sd"], sqrt(Fit$xiError[, 1]), xlab = "RMSE (A)", 
            ylab = "IRT Standard Error", pch = 19, main = "A: IRT SE vs. RMSE")
        abline(a = 0, b = 1)
        plot(RMSE$stats[RMSE$stats$N == examinees & RMSE$stats$J == 
            items, "B_sd"], sqrt(Fit$xiError[, 2]), xlab = "RMSE (B)", 
            ylab = "IRT Standard Error", pch = 19, main = "B: IRT SE vs. RMSE")
        abline(a = 0, b = 1)
    }
    else {
        par(mfrow = c(2, 3))
        plot(Fit$xi[, 1], Fit$EmpSE$MCSA[, 1], xlab = "A1 (Generated)", 
            ylab = "Standard Error", pch = 19, main = "dA1: Comparing initseq() SE to RMSE")
        abline(lm(Fit$EmpSE$MCSA[, 1] ~ Fit$xi[, 1])$coefficients, 
            lwd = 1.5)
        points(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, c("A1", "A1_sd")], col = 2, pch = 19)
        abline(lm(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, c("A1_sd", "A1")])$coefficients, lwd = 1.5, 
            col = 2)
        legend("topleft", c("MCMC", "RMSE"), pch = 19, col = c(1, 
            2))
        plot(Fit$xi[, 2], Fit$EmpSE$MCSA[, 2], xlab = "A2 (Generated)", 
            ylab = "Standard Error", pch = 19, main = "dA2: Comparing initseq() SE to RMSE")
        abline(lm(Fit$EmpSE$MCSA[, 2] ~ Fit$xi[, 2])$coefficients, 
            lwd = 1.5)
        points(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, c("A2", "A2_sd")], col = 2, pch = 19)
        abline(lm(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, c("A2_sd", "A2")])$coefficients, lwd = 1.5, 
            col = 2)
        legend("topleft", c("MCMC", "RMSE"), pch = 19, col = c(1, 
            2))
        plot(Fit$B, Fit$EmpSE$MCSB, xlab = "B (Generated)", ylab = "Standard Error", 
            pch = 19, main = "dB: Comparing initseq() SE to RMSE")
        lines(Fit$B[order(Fit$B)], predict(lm(Fit$EmpSE$MCSB ~ 
            Fit$B + I(Fit$B^2)))[order(Fit$B)], lwd = 1.5)
        points(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, c("B", "B_sd")], col = 2, pch = 19)
        lines(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, "B"][order(RMSE$STATS[RMSE$STATS$N == examinees & 
            RMSE$STATS$J == items, "B"])], predict(lm(RMSE$STATS[RMSE$STATS$N == 
            examinees & RMSE$STATS$J == items, "B_sd"] ~ RMSE$STATS[RMSE$STATS$N == 
            examinees & RMSE$STATS$J == items, "B"] + I(RMSE$STATS[RMSE$STATS$N == 
            examinees & RMSE$STATS$J == items, "B"]^2)))[order(RMSE$STATS[RMSE$STATS$N == 
            examinees & RMSE$STATS$J == items, "B"])], lwd = 1.5, 
            col = 2)
        legend("topleft", c("MCMC", "RMSE"), pch = 19, col = c(1, 
            2))
        plot(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, "A1_sd"], Fit$EmpSE$MCSA[, 1], xlab = "RMSE (A1)", 
            ylab = "MCMC Standard Error", pch = 19, main = "A1: MCMC SE vs. RMSE")
        abline(a = 0, b = 1)
        plot(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, "A2_sd"], Fit$EmpSE$MCSA[, 2], xlab = "RMSE (A2)", 
            ylab = "MCMC Standard Error", pch = 19, main = "A2: MCMC SE vs. RMSE")
        abline(a = 0, b = 1)
        plot(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, "B_sd"], Fit$EmpSE$MCSB, xlab = "RMSE (B)", 
            ylab = "MCMC Standard Error", pch = 19, main = "B: MCMC SE vs. RMSE")
        abline(a = 0, b = 1)
        par(mfrow = c(2, 3))
        plot(Fit$xi[, 1], sqrt(Fit$EmpSE$VARLMIo[, 1]), xlab = "A1 (Generated)", 
            ylab = "Standard Error", pch = 19, main = "dA1: Comparing Ogive Emp SE to RMSE")
        abline(lm(sqrt(Fit$EmpSE$VARLMIo[, 1]) ~ Fit$xi[, 1])$coefficients, 
            lwd = 1.5)
        points(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, c("A1", "A1_sd")], col = 2, pch = 19)
        abline(lm(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, c("A1_sd", "A1")])$coefficients, lwd = 1.5, 
            col = 2)
        legend("topleft", c("Ogive Emp SE", "RMSE"), pch = 19, 
            col = c(1, 2))
        plot(Fit$xi[, 2], sqrt(Fit$EmpSE$VARLMIo[, 2]), xlab = "A2 (Generated)", 
            ylab = "Standard Error", pch = 19, main = "dA2: Comparing Ogive Emp SE to RMSE")
        abline(lm(sqrt(Fit$EmpSE$VARLMIo[, 2]) ~ Fit$xi[, 2])$coefficients, 
            lwd = 1.5)
        points(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, c("A2", "A2_sd")], col = 2, pch = 19)
        abline(lm(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, c("A2_sd", "A2")])$coefficients, lwd = 1.5, 
            col = 2)
        legend("topleft", c("Ogive Emp SE", "RMSE"), pch = 19, 
            col = c(1, 2))
        plot(Fit$B, sqrt(abs(Fit$EmpSE$VARLMIo[, 3])), xlab = "B (Generated)", 
            ylab = "Standard Error", pch = 19, main = "dB: Comparing Ogive Emp SE to RMSE")
        lines(Fit$B[order(Fit$B)], predict(lm(sqrt(abs(Fit$EmpSE$VARLMIo[, 
            3])) ~ Fit$B + I(Fit$B^2)))[order(Fit$B)], lwd = 1.5)
        points(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, c("B", "B_sd")], col = 2, pch = 19)
        lines(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, "B"][order(RMSE$STATS[RMSE$STATS$N == examinees & 
            RMSE$STATS$J == items, "B"])], predict(lm(RMSE$STATS[RMSE$STATS$N == 
            examinees & RMSE$STATS$J == items, "B_sd"] ~ RMSE$STATS[RMSE$STATS$N == 
            examinees & RMSE$STATS$J == items, "B"] + I(RMSE$STATS[RMSE$STATS$N == 
            examinees & RMSE$STATS$J == items, "B"]^2)))[order(RMSE$STATS[RMSE$STATS$N == 
            examinees & RMSE$STATS$J == items, "B"])], lwd = 1.5, 
            col = 2)
        legend("topleft", c("Ogive Emp SE", "RMSE"), pch = 19, 
            col = c(1, 2))
        plot(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, "A1_sd"], sqrt(Fit$EmpSE$VARLMIo[, 1]), xlab = "RMSE (A1)", 
            ylab = "Ogive Emp Standard Error", pch = 19, main = "A1: Ogive Emp SE vs. RMSE")
        abline(a = 0, b = 1)
        plot(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, "A2_sd"], sqrt(Fit$EmpSE$VARLMIo[, 2]), xlab = "RMSE (A2)", 
            ylab = "Ogive Emp Standard Error", pch = 19, main = "A2: Ogive Emp SE vs. RMSE")
        abline(a = 0, b = 1)
        plot(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, "B_sd"], sqrt(Fit$EmpSE$VARLMIo[, 3]), xlab = "RMSE (B)", 
            ylab = "Ogive Emp Standard Error", pch = 19, main = "B: Ogive Emp SE vs. RMSE")
        abline(a = 0, b = 1)
        par(mfrow = c(2, 3))
        plot(Fit$xi[, 1], sqrt(Fit$EmpSE$VARLMIi[, 1]), xlab = "A1 (Generated)", 
            ylab = "Standard Error", pch = 19, main = "dA1: Comparing IRT Emp SE to RMSE")
        abline(lm(sqrt(Fit$EmpSE$VARLMIi[, 1]) ~ Fit$xi[, 1])$coefficients, 
            lwd = 1.5)
        points(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, c("A1", "A1_sd")], col = 2, pch = 19)
        abline(lm(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, c("A1_sd", "A1")])$coefficients, lwd = 1.5, 
            col = 2)
        legend("topleft", c("IRT Emp SE", "RMSE"), pch = 19, 
            col = c(1, 2))
        plot(Fit$xi[, 2], sqrt(Fit$EmpSE$VARLMIi[, 2]), xlab = "A2 (Generated)", 
            ylab = "Standard Error", pch = 19, main = "dA2: Comparing IRT Emp SE to RMSE")
        abline(lm(sqrt(Fit$EmpSE$VARLMIi[, 2]) ~ Fit$xi[, 2])$coefficients, 
            lwd = 1.5)
        points(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, c("A2", "A2_sd")], col = 2, pch = 19)
        abline(lm(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, c("A2_sd", "A2")])$coefficients, lwd = 1.5, 
            col = 2)
        legend("topleft", c("IRT Emp SE", "RMSE"), pch = 19, 
            col = c(1, 2))
        plot(Fit$B, sqrt(abs(Fit$EmpSE$VARLMIi[, 3])), xlab = "B (Generated)", 
            ylab = "Standard Error", pch = 19, main = "dB: Comparing IRT Emp SE to RMSE")
        lines(Fit$B[order(Fit$B)], predict(lm(sqrt(abs(Fit$EmpSE$VARLMIi[, 
            3])) ~ Fit$B + I(Fit$B^2)))[order(Fit$B)], lwd = 1.5)
        points(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, c("B", "B_sd")], col = 2, pch = 19)
        lines(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, "B"][order(RMSE$STATS[RMSE$STATS$N == examinees & 
            RMSE$STATS$J == items, "B"])], predict(lm(RMSE$STATS[RMSE$STATS$N == 
            examinees & RMSE$STATS$J == items, "B_sd"] ~ RMSE$STATS[RMSE$STATS$N == 
            examinees & RMSE$STATS$J == items, "B"] + I(RMSE$STATS[RMSE$STATS$N == 
            examinees & RMSE$STATS$J == items, "B"]^2)))[order(RMSE$STATS[RMSE$STATS$N == 
            examinees & RMSE$STATS$J == items, "B"])], lwd = 1.5, 
            col = 2)
        legend("topleft", c("IRT Emp SE", "RMSE"), pch = 19, 
            col = c(1, 2))
        plot(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, "A1_sd"], sqrt(Fit$EmpSE$VARLMIi[, 1]), xlab = "RMSE (A1)", 
            ylab = "IRT Emp Standard Error", pch = 19, main = "A1: IRT Emp SE vs. RMSE")
        abline(a = 0, b = 1)
        plot(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, "A2_sd"], sqrt(Fit$EmpSE$VARLMIi[, 2]), xlab = "RMSE (A2)", 
            ylab = "IRT Emp Standard Error", pch = 19, main = "A2: IRT Emp SE vs. RMSE")
        abline(a = 0, b = 1)
        plot(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, "B_sd"], sqrt(Fit$EmpSE$VARLMIi[, 3]), xlab = "RMSE (B)", 
            ylab = "IRT Emp Standard Error", pch = 19, main = "B: IRT Emp SE vs. RMSE")
        abline(a = 0, b = 1)
        par(mfrow = c(2, 3))
        plot(Fit$xi[, 1], sqrt(abs(Fit$xiError[, 1])), xlab = "A1 (Generated)", 
            ylab = "Standard Error", pch = 19, main = "dA1: Comparing IRT SE to RMSE")
        abline(lm(sqrt(abs(Fit$xiError[, 1])) ~ Fit$xi[, 1])$coefficients, 
            lwd = 1.5)
        points(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, c("A1", "A1_sd")], col = 2, pch = 19)
        abline(lm(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, c("A1_sd", "A1")])$coefficients, lwd = 1.5, 
            col = 2)
        legend("topleft", c("IRT SE", "RMSE"), pch = 19, col = c(1, 
            2))
        plot(Fit$xi[, 2], sqrt(abs(Fit$xiError[, 2])), xlab = "A2 (Generated)", 
            ylab = "Standard Error", pch = 19, main = "dA2: Comparing IRT SE to RMSE")
        abline(lm(sqrt(abs(Fit$xiError[, 2])) ~ Fit$xi[, 2])$coefficients, 
            lwd = 1.5)
        points(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, c("A2", "A2_sd")], col = 2, pch = 19)
        abline(lm(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, c("A2_sd", "A2")])$coefficients, lwd = 1.5, 
            col = 2)
        legend("topleft", c("IRT SE", "RMSE"), pch = 19, col = c(1, 
            2))
        plot(Fit$B, sqrt(abs(Fit$xiError[, 3])), xlab = "B (Generated)", 
            ylab = "Standard Error", pch = 19, main = "dB: Comparing IRT SE to RMSE")
        lines(Fit$B[order(Fit$B)], predict(lm(sqrt(abs(Fit$xiError[, 
            3])) ~ Fit$B + I(Fit$B^2)))[order(Fit$B)], lwd = 1.5)
        points(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, c("B", "B_sd")], col = 2, pch = 19)
        lines(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, "B"][order(RMSE$STATS[RMSE$STATS$N == examinees & 
            RMSE$STATS$J == items, "B"])], predict(lm(RMSE$STATS[RMSE$STATS$N == 
            examinees & RMSE$STATS$J == items, "B_sd"] ~ RMSE$STATS[RMSE$STATS$N == 
            examinees & RMSE$STATS$J == items, "B"] + I(RMSE$STATS[RMSE$STATS$N == 
            examinees & RMSE$STATS$J == items, "B"]^2)))[order(RMSE$STATS[RMSE$STATS$N == 
            examinees & RMSE$STATS$J == items, "B"])], lwd = 1.5, 
            col = 2)
        legend("topleft", c("IRT SE", "RMSE"), pch = 19, col = c(1, 
            2))
        plot(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, "A1_sd"], sqrt(Fit$xiError[, 1]), xlab = "RMSE (A1)", 
            ylab = "IRT Standard Error", pch = 19, main = "A1: IRT SE vs. RMSE")
        abline(a = 0, b = 1)
        plot(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, "A2_sd"], sqrt(Fit$xiError[, 2]), xlab = "RMSE (A2)", 
            ylab = "IRT Standard Error", pch = 19, main = "A2: IRT SE vs. RMSE")
        abline(a = 0, b = 1)
        plot(RMSE$STATS[RMSE$STATS$N == examinees & RMSE$STATS$J == 
            items, "B_sd"], sqrt(Fit$xiError[, 3]), xlab = "RMSE (B)", 
            ylab = "IRT Standard Error", pch = 19, main = "B: IRT SE vs. RMSE")
        abline(a = 0, b = 1)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
