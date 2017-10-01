\name{Info3PL}
\alias{Info3PL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Info3PL(aa, bb, cc, Post)
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
  \item{Post}{
%%     ~~Describe \code{Post} here~~
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
function (aa, bb, cc, Post) 
{
    NJ <- mat.or.vec(length(Post), 1) + 1
    Z <- outer(Post, bb, "-")
    print("dim Z")
    print(dim(Z))
    P <- threePL(aa, bb, cc, Post)
    print("dim P")
    print(dim(P))
    A <- t(aa \%*\% t(NJ))
    print("dim A")
    print(dim(A))
    C <- t(cc \%*\% t(NJ))
    print("dim C")
    print(dim(C))
    IA <- (1.7/(1 - C))^2 * Z^2 * (P - C)^2 * (1 - P)/P
    IB <- ((1.7 * A)/(1 - C))^2 * (P - C)^2 * (1 - P)/P
    IC <- 1/(1 - C)^2 * (1 - P)/P
    IAB <- (-1) * IA * A/Z
    IAC <- IC * 1.7 * Z * (P - C)
    IBC <- IC * (-1.7) * A * (P - C)
    IA <- colSums(IA)
    IB <- colSums(IB)
    IC <- colSums(IC)
    IAB <- colSums(IAB)
    IAC <- colSums(IAC)
    IBC <- colSums(IBC)
    IM <- matrix(c(IA, IB, IC, IAB, IAC, IBC), nrow = length(aa), 
        ncol = 6)
    print("Information Matrix Done")
    print(IM)
    return(IM)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
