\name{CheckParams}
\alias{CheckParams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check Parameters for Generation or Fit
}
\description{
Verifies names of parameters, cleans case, creates defaults if parameters not defined.
}
\usage{
CheckParams(parameters,generate=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parameters}{
list; the list of parameters for generating or fitting
}
  \item{generate}{
TRUE/FALSE; if true, parameters are assumed to be the "generate" parameters, else, the model's parameters for fitting.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a parameter list
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

structure<-CheckParams()
settings<-CheckParams(generate=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
