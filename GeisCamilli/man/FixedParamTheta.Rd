\name{FixedParamTheta}
\alias{FixedParamTheta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FixedParamTheta(FitDATA, rp, IT = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FitDATA}{
%%     ~~Describe \code{FitDATA} here~~
}
  \item{rp}{
%%     ~~Describe \code{rp} here~~
}
  \item{IT}{
%%     ~~Describe \code{IT} here~~
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
function (FitDATA, rp, IT = 100) 
{
    if (settings$Adim > 1) {
        if (!is.na(FitDATA$AR)) {
            A = as.matrix(FitDATA$AR$loadings)
            THat = as.matrix(FitDATA$Trot[, 1:settings$Adim])
        }
        else {
            A = as.matrix(FitDATA$A)
            THat = as.matrix(FitDATA$That[, 1:settings$Adim])
        }
        B = FitDATA$B
        C = FitDATA$C
    }
    else {
        A = FitDATA$A
        B = FitDATA$B
        C = FitDATA$C
        THat = FitDATA$That[, 1]
    }
    Alines <- A
    Blines <- B
    settings = FitDATA$settings
    if (settings$fm \%in\% c("new", "pca", "eigen")) 
        settings$fm <- "camilli"
    prior <- list(tmu = settings$tmu, tsigma = settings$tsigma)
    J <- ncol(rp)
    N <- nrow(rp)
    ifelse(settings$guess, W <- DrawW(aa = A, bb = B, cc = C, 
        tt = THat, rp = rp), W <- NA)
    Z <- SampZ(aa = A, bb = B, that = THat, rp = rp, w = W)
    AFSEiter <- array(A, dim = c(J, settings$Adim, 1))
    BFSEiter <- matrix(B, nrow = J, ncol = 1)
    ifelse(settings$guess, CFSEiter <- matrix(C, nrow = J, ncol = 1), 
        CFSEiter <- NA)
    FTiter <- array(THat, dim = c(N, settings$Adim, 1))
    for (it in 1:IT) {
        Z <- SampZ(aa = A, bb = B, that = THat, rp = rp, w = W)
        LL <- GIFAFullLL(A, B, Z, THat, prior = prior)
        FTiter <- abind(FTiter, as.matrix(THat), along = 3)
        if (settings$guess) {
            PSI <- GIFAEstimate(aa = A, bb = B, zz = Z, tt = THat, 
                settings = settings, w = W, rp = rp)
        }
        else {
            PSI <- GIFAEstimate(aa = A, bb = B, zz = Z, tt = THat, 
                settings = settings)
        }
        THat <- SampT(aa = PSI$A, bb = PSI$B, zz = Z, rp = rp, 
            prior = prior)
        AFSEiter <- abind(AFSEiter, as.matrix(PSI$A), along = 3)
        BFSEiter <- cbind(BFSEiter, PSI$B)
        if (settings$guess) 
            CFSEiter <- cbind(CFSEiter, PSI$C)
    }
    if (settings$Adim > 1) {
        FTErr <- as.matrix(apply(FTiter, c(1, 2), sd))
        FAErr <- as.matrix(apply(AFSEiter, c(1, 2), sd))
        FT <- as.matrix(apply(FTiter, c(1, 2), mean))
    }
    else {
        FTErr <- as.vector(apply(as.matrix(FTiter[, 1, ]), 1, 
            sd))
        FAErr <- as.matrix(apply(AFSEiter[, 1, ], 1, sd))
        FT <- as.vector(apply(as.matrix(FTiter[, 1, ]), 1, mean))
    }
    FBErr <- as.vector(apply(BFSEiter, 1, sd))
    if (settings$plots) {
        if (settings$Adim > 1) {
            par(mfrow = c(1, settings$Adim), mar = c(3, 3, 1, 
                1))
            for (i in 1:settings$Adim) {
                plot(FTiter[1, i, ], type = "n", xlim = c(0, 
                  length(FTiter[1, 1, ])), ylim = range(FTiter), 
                  main = paste("Theta", i, "of 10 Examinees"))
                for (j in 1:10) {
                  lines(FTiter[j, i, ], col = j)
                }
                abline(h = FT[1:10, i])
            }
            par(mfrow = c(1, settings$Adim + 1), mar = c(3, 3, 
                1, 1))
            for (i in 1:settings$Adim) {
                plot(AFSEiter[1, i, ], type = "n", xlim = c(0, 
                  length(AFSEiter[1, 1, ])), ylim = range(AFSEiter), 
                  main = paste("A", i))
                for (j in 1:J) {
                  lines(AFSEiter[j, i, ], col = j)
                }
                abline(h = Alines[, i])
            }
            plot(BFSEiter[1, ], type = "n", xlim = c(0, length(BFSEiter[1, 
                ])), ylim = range(BFSEiter), main = "B")
            for (j in 1:J) {
                lines(BFSEiter[j, ], col = j)
            }
            abline(h = Blines)
        }
        else {
            par(mfrow = c(1, 1), mar = c(3, 3, 1, 1))
            plot(FTiter[1, 1, ], type = "n", xlim = c(0, length(FTiter[1, 
                1, ])), ylim = range(FTiter), main = paste("Simulating Theta of 10 Examinees"))
            for (j in 1:10) {
                lines(FTiter[j, 1, ], col = j)
            }
            abline(h = FT[1:10])
            par(mfrow = c(1, 2), mar = c(3, 3, 1, 1))
            plot(AFSEiter[1, 1, ], type = "n", xlim = c(0, length(AFSEiter[1, 
                1, ])), ylim = range(AFSEiter), main = paste("A", 
                1))
            for (j in 1:J) {
                lines(AFSEiter[j, 1, ], col = j)
            }
            abline(h = Alines)
            plot(BFSEiter[1, ], type = "n", xlim = c(0, length(BFSEiter[1, 
                ])), ylim = range(BFSEiter), main = "B")
            for (j in 1:J) {
                lines(BFSEiter[j, ], col = j)
            }
            abline(h = Blines)
        }
    }
    return(list(FT = FT, FTErr = FTErr, FAErr = FAErr, FBErr = FBErr))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
