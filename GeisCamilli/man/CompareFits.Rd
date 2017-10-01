\name{CompareFits}
\alias{CompareFits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CompareFits(Gen = NA, Fit1, Fit2, Fit3 = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Gen}{
%%     ~~Describe \code{Gen} here~~
}
  \item{Fit1}{
%%     ~~Describe \code{Fit1} here~~
}
  \item{Fit2}{
%%     ~~Describe \code{Fit2} here~~
}
  \item{Fit3}{
%%     ~~Describe \code{Fit3} here~~
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
function (Gen = NA, Fit1, Fit2, Fit3 = NA) 
{
    Err1 <- abs(Fit1$xiError)^0.5
    Err2 <- abs(Fit2$xiError)^0.5
    if (!is.na(Fit3)) 
        Err3 <- abs(Fit3$xiError)^0.5
    UBnds <- apply(rbind(Fit1$xi, Fit2$xi), 2, max)
    LBnds <- apply(rbind(Fit1$xi, Fit2$xi), 2, min)
    plot(Fit1$xi[, 1], Fit1$xi[, 2], pch = 1, cex = 1.2, col = "red", 
        main = "B vs. A", xlim = c(LBnds[1] - 0.2, UBnds[1] + 
            0.2), ylim = c(LBnds[2] - 0.2, UBnds[2] + 0.2))
    arrows(Fit1$xi[, 1] - 2 * Err1[, 1], Fit1$xi[, 2], Fit1$xi[, 
        1] + 2 * Err1[, 1], Fit1$xi[, 2], code = 3, col = "red", 
        angle = 90, length = 0.07)
    arrows(Fit1$xi[, 1], Fit1$xi[, 2] - 2 * Err1[, 2], Fit1$xi[, 
        1], Fit1$xi[, 2] + 2 * Err1[, 2], code = 3, col = "red", 
        angle = 90, length = 0.07)
    points(Fit2$xi[, 1], Fit2$xi[, 2], pch = 19, col = "blue")
    arrows(Fit2$xi[, 1] - 2 * Err2[, 1], Fit2$xi[, 2], Fit2$xi[, 
        1] + 2 * Err2[, 1], Fit2$xi[, 2], code = 3, col = "blue", 
        angle = 90, length = 0.07)
    arrows(Fit2$xi[, 1], Fit2$xi[, 2] - 2 * Err2[, 2], Fit2$xi[, 
        1], Fit2$xi[, 2] + 2 * Err2[, 2], code = 3, col = "blue", 
        angle = 90, length = 0.07)
    legend("topright", c("SAEM", "EM"), pch = c(1, 19), col = c("red", 
        "blue"))
    if (!is.na(Gen)) {
        points(Gen$XI[, 1], Gen$XI[, 2], pch = 13)
        AX1 <- Fit1$xi[, 1] - Gen$XI[, 1]
        BX1 <- Fit1$xi[, 2] - Gen$XI[, 2]
        AX2 <- Fit2$xi[, 1] - Gen$XI[, 1]
        BX2 <- Fit2$xi[, 2] - Gen$XI[, 2]
        plot(c(AX1, BX1), pch = 1, cex = 1.2, col = "red", main = "Errors", 
            ylim = 1.5 * range(c(AX1, BX1)))
        arrows(1:length(AX1), AX1 - 2 * Err1[, 1], 1:length(AX1), 
            AX1 + 2 * Err1[, 1], code = 3, col = "red", angle = 90, 
            length = 0.07)
        arrows(length(AX1) + 1:length(BX1), BX1 - 2 * Err1[, 
            2], length(AX1) + 1:length(BX1), BX1 + 2 * Err1[, 
            2], code = 3, col = "red", angle = 90, length = 0.07)
        points(1:length(c(AX1, BX1)), c(AX2, BX2), pch = 19, 
            col = "blue")
        arrows(1:length(AX1), AX2 - 2 * Err2[, 1], 1:length(AX1), 
            AX2 + 2 * Err2[, 1], code = 3, col = "blue", angle = 90, 
            length = 0.07)
        arrows(length(AX1) + 1:length(BX1), BX2 - 2 * Err2[, 
            2], length(AX1) + 1:length(BX1), BX2 + 2 * Err2[, 
            2], code = 3, col = "blue", angle = 90, length = 0.07)
        legend("topright", c("SAEM", "EM"), pch = c(1, 19), col = c("red", 
            "blue"))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
