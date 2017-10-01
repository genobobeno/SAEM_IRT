\name{SaveRunInfo}
\alias{SaveRunInfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SaveRunInfo(gen.xi, gen.theta, gen.RP, xi, xiErr, theta, thetaErr, settings = settings, structure = structure, tag)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gen.xi}{
%%     ~~Describe \code{gen.xi} here~~
}
  \item{gen.theta}{
%%     ~~Describe \code{gen.theta} here~~
}
  \item{gen.RP}{
%%     ~~Describe \code{gen.RP} here~~
}
  \item{xi}{
%%     ~~Describe \code{xi} here~~
}
  \item{xiErr}{
%%     ~~Describe \code{xiErr} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{thetaErr}{
%%     ~~Describe \code{thetaErr} here~~
}
  \item{settings}{
%%     ~~Describe \code{settings} here~~
}
  \item{structure}{
%%     ~~Describe \code{structure} here~~
}
  \item{tag}{
%%     ~~Describe \code{tag} here~~
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
function (gen.xi, gen.theta, gen.RP, xi, xiErr, theta, thetaErr, 
    settings = settings, structure = structure, tag) 
{
    SimData = list(gen.xi = gen.xi, gen.theta = gen.theta, gen.RP = gen.RP, 
        xi = xi, xiErr = xiErr, theta = theta, thetaErr = thetaErr, 
        settings = settings, structure = structure)
    save(SimData, file = paste("D", settings$dimF, "_J", nrow(xi), 
        "_N", nrow(as.matrix(theta)), "_", tag, ".Rda", sep = ""))
    print(paste("Saved file: D", settings$dimF, "_J", nrow(xi), 
        "_N", nrow(as.matrix(theta)), "_", tag, ".Rda", sep = ""))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
