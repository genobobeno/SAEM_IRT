\name{GetThetaHat}
\alias{GetThetaHat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GetThetaHat(aa, bb, cc, rp, tHat, zHat, w, prior, setting, R = NA)
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
  \item{tHat}{
%%     ~~Describe \code{tHat} here~~
}
  \item{zHat}{
%%     ~~Describe \code{zHat} here~~
}
  \item{w}{
%%     ~~Describe \code{w} here~~
}
  \item{prior}{
%%     ~~Describe \code{prior} here~~
}
  \item{setting}{
%%     ~~Describe \code{setting} here~~
}
  \item{R}{
%%     ~~Describe \code{R} here~~
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
function (aa, bb, cc, rp, tHat, zHat, w, prior, setting, R = NA) 
{
    ZHAT <- zHat
    if (setting$Adim > 1) {
        THAT <- array(tHat, dim = c(nrow(rp), setting$Adim, 1))
    }
    else {
        THAT <- tHat
    }
    for (i in 1:setting$nesttheta) {
        tHat <- SampT(aa = aa, bb = bb, zz = zHat, rp = rp, prior = prior)
        zHat <- SampZ(aa = aa, bb = bb, that = tHat, rp = rp, 
            w = w)
        ZHAT <- cbind(ZHAT, zHat)
        ifelse(setting$Adim > 1, THAT <- abind(THAT, tHat, along = 3), 
            THAT <- cbind(THAT, tHat))
    }
    tHat <- SampT(aa = aa, bb = bb, zz = zHat, rp = rp, prior = prior)
    ifelse(setting$Adim > 1, THAT <- abind(THAT[, , -1], tHat, 
        along = 3), THAT <- cbind(THAT[, -1], tHat))
    if (setting$Adim > 1) {
        THETA <- apply(THAT, c(1, 2), mean)
        THETA <- cbind(THETA, THETA - apply(THAT, c(1, 2), sd), 
            THETA + apply(THAT, c(1, 2), sd))
        colnames(THETA) <- paste(c("Theta", "Theta-SE", "Theta+SE"), 
            rep(c(1:setting$Adim), rep(3, setting$Adim)))
    }
    else {
        THETA <- rowMeans(THAT)
        THETA <- cbind(THETA, THETA - apply(THAT, 1, sd), THETA + 
            apply(THAT, 1, sd))
        colnames(THETA) <- c("Theta", "Theta-SE", "Theta+SE")
    }
    if (setting$thetamap) {
        TMAP <- ThetaMAP(aa = aa, bb = bb, cc = cc, rp = rp, 
            settings = setting)
        TMAP <- as.matrix(TMAP)
        ifelse(setting$Adim == 1, colnames(TMAP) <- "TMAP", colnames(TMAP) <- paste("TMAP", 
            1:setting$Adim, sep = ""))
        if (!is.na(R)) {
            aa = R$loadings
            TRMAP <- ThetaMAP(aa = aa, bb = bb, cc = cc, rp = rp, 
                settings = setting)
            TRMAP <- as.matrix(TRMAP)
            ifelse(setting$Adim == 1, colnames(TRMAP) <- "TRMAP", 
                colnames(TRMAP) <- paste("TRMAP", 1:setting$Adim, 
                  sep = ""))
        }
        else {
            TRMAP = NA
        }
    }
    else {
        TMAP <- NA
        TRMAP = NA
    }
    return(list(THETA = THETA, TMAP = TMAP, TRMAP = TRMAP))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
