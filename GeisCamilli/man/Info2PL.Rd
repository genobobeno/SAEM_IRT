\name{Info2PL}
\alias{Info2PL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Info2PL(aa, bb, Post)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aa}{
%%     ~~Describe \code{aa} here~~
}
  \item{bb}{
%%     ~~Describe \code{bb} here~~
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
function (aa, bb, Post) 
{
    NJ <- mat.or.vec(length(Post), 1) + 1
    Z <- outer(Post, bb, "-")
    print("dim Z")
    print(dim(Z))
    P <- twoPL(aa, bb, Post)
    print("dim P")
    print(dim(P))
    A <- t(aa \%*\% t(NJ))
    print("dim A")
    print(dim(A))
    IA <- (1.7)^2 * Z^2 * (P) * (1 - P)
    IB <- (1.7 * A)^2 * (P) * (1 - P)
    IAB <- (-1) * IA * A/Z
    IA <- colSums(IA)
    IB <- colSums(IB)
    IAB <- colSums(IAB)
    IM <- matrix(c(IA, IB, IAB), nrow = length(aa), ncol = 3)
    print("Information Matrix Done")
    print(IM)
    return(IM)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
