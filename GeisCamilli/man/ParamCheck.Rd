\name{ParamCheck}
\alias{ParamCheck}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check Parameters for Generation or Fit
}
\description{
Verifies names of parameters, cleans case, creates defaults if parameters not defined.
}
\usage{
ParamCheck(name, plist, dplist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{
Name of default parameter
}
  \item{plist}{
Parameter list defined by user
}
  \item{dplist}{
Default parameter list defined by package; 'structure' for Generator, 'settings' for Fit.
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

## The function is currently defined as
function (name, plist, dplist) 
{
    if (!(tolower(name) \%in\% tolower(names(dplist)) && !is.na(plist[name]))) {
        plist[name] = dplist[name]
    }
    names(plist)[names(plist) == name] <- names(dplist)[which(tolower(names(dplist)) \%in\% 
        tolower(name))]
    return(plist)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
