\name{Timing}
\alias{Timing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Start Time or Calculate Time Since Start
}
\description{
Starts a timer, or tells you how long it's been since you started the timer.
}
\usage{
Timing(start=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{start}{If NA, starts timer. If a time, tells you duration.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
returns time.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Eugene Geis
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
function (settings) 
{
    if (tolower(settings$quad == "manual")) {
        Grid <- seq(settings$gridbounds[1], settings$gridbounds[2], 
            length.out = settings$nodes)
        return(cbind(Grid, dnorm(Grid, mean = settings$tmu, sd = settings$tsigma)))
    }
    else if (tolower(settings$quad == "gauss-hermite")) {
        Grid <- gaussHermiteData(settings$nodes)
        return(cbind(Grid$x, Grid$w))
    }
    else {
        print("Quadrature method not chosen, or not among choices")
        print("Please choose manual, or gauss-hermite")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
