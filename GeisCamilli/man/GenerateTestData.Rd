\name{GenerateTestData}
\alias{GenerateTestData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GenerateTestData(j = 30, n = 1000, structure = structure)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{j}{
%%     ~~Describe \code{j} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{structure}{
%%     ~~Describe \code{structure} here~~
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
function (j = 30, n = 1000, structure = structure) 
{
    stopifnot(tolower(structure$icc) \%in\% c("logistic", "ogive"), 
        j > 4, length(structure$Aparams) >= 1, n > 50, structure$Adim == 
            length(structure$tmu))
    print("***********************  Generating Parameters  **************************")
    A <- GenerateA(j, structure$Adim, tolower(structure$Adist), 
        structure$Aparams)
    b <- GenerateB(j, tolower(structure$bdist), structure$bparams)
    if (length(A) == j) 
        cat(" ...IRT pseudo-difficulties of Ogive: ", b/A)
    xi <- cbind(A, b)
    if (structure$guess) {
        c <- GenerateC(j, tolower(structure$cdist), structure$cparams)
        xi <- cbind(xi, c)
    }
    t <- GenerateTheta(n, structure$tmu, structure$tsigma)
    rp <- GenerateRP(xi, t, structure = structure)
    gen.rp <- rp
    gen.xi <- xi
    gen.theta <- t
    gen.structure <- structure
    if (!is.na(structure$simfile)) {
        if (grepl("\\.[Rr][Dd][Aa]", structure$simfile)) {
            save(gen.rp, gen.xi, gen.theta, gen.structure, file = structure$simfile)
        }
        else {
            save(gen.rp, gen.xi, gen.theta, gen.structure, file = paste(structure$simfile, 
                ".rda", sep = ""))
        }
    }
    return(list(RP = rp, XI = xi, THETA = t))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
