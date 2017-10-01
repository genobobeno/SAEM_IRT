\name{GenerateA}
\alias{GenerateA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GenerateA(j, Adim, Adist, Aparams)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{j}{
%%     ~~Describe \code{j} here~~
}
  \item{Adim}{
%%     ~~Describe \code{Adim} here~~
}
  \item{Adist}{
%%     ~~Describe \code{Adist} here~~
}
  \item{Aparams}{
%%     ~~Describe \code{Aparams} here~~
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
function (j, Adim, Adist, Aparams) 
{
    if (Adist == "unif") {
        stopifnot(length(Aparams) == 2)
        A <- runif(j, min = Aparams[1], max = Aparams[2])
        if (Adim > 1) {
            for (i in 2:Adim) {
                A <- cbind(A, runif(j, min = Aparams[1] - 0.1, 
                  max = 1))
                if (i == 2) 
                  A[1:((Adim - 1) * ceiling(j/Adim)), i] <- 0
                if (i >= 3) {
                  A[1:((Adim - (i - 1)) * ceiling(j/Adim)), i] <- 0
                  A[(((Adim - (i - 2)) * ceiling(j/Adim)) + 1):j, 
                    i] <- 0
                }
            }
        }
    }
    else if (Adist == "norm") {
        stopifnot(length(Aparams) == 2, Aparams[1] > 0.2, Aparams[2] < 
            1)
        A <- rtruncnorm(j, a = 0.2, b = 2, Aparams[1], Aparams[2])
        if (Adim > 1) {
            for (i in 2:Adim) {
                A <- cbind(A, rtruncnorm(j, a = 0.2, b = 2, Aparams[1], 
                  Aparams[2]))
                if (i == 2) 
                  A[1:((Adim - 1) * ceiling(j/Adim)), i] <- 0
                if (i >= 3) {
                  A[1:((Adim - (i - 1)) * ceiling(j/Adim)), i] <- 0
                  A[(((Adim - (i - 2)) * ceiling(j/Adim)) + 1):j, 
                    i] <- 0
                }
            }
        }
    }
    else if (Adist == "beta") {
        stopifnot(length(Aparams) == 2, Aparams[1] > 0, Aparams[2] < 
            2.5)
        A <- rbeta(j, 2.5, 3)
        A <- (Aparams[2] - Aparams[1]) * A + Aparams[1]
        if (Adim > 1) {
            for (i in 2:Adim) {
                A <- cbind(A, (1 - 0.1) * rbeta(j, 2.5, 3) + 
                  0.1)
                if (i == 2) 
                  A[1:((Adim - 1) * ceiling(j/Adim)), i] <- 0
                if (i >= 3) {
                  A[1:((Adim - (i - 1)) * ceiling(j/Adim)), i] <- 0
                  A[(((Adim - (i - 2)) * ceiling(j/Adim)) + 1):j, 
                    i] <- 0
                }
            }
        }
    }
    return(A)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
