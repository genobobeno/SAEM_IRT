\name{ReadData}
\alias{ReadData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ReadData(GParam = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GParam}{
%%     ~~Describe \code{GParam} here~~
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
function (GParam = FALSE) 
{
    TTrue <- as.matrix(read.table("Thetas.dat", header = F))
    A_Gen <- as.matrix(read.table("Aparams.dat", header = F))
    B_Gen <- as.matrix(read.table("Bparams.dat", header = F))
    if (GParam == TRUE) {
        C_Gen <- as.matrix(read.table("Cparams.dat", header = F))
    }
    Response <- read.table("Response.dat", header = F)
    if (GParam == FALSE) {
        RList <- list(A_Gen = A_Gen, B_Gen = B_Gen, Response = Response, 
            TTrue = TTrue)
    }
    else {
        RList <- list(A_Gen = A_Gen, B_Gen = B_Gen, C_Gen = C_Gen, 
            Response = Response, TTrue = TTrue)
    }
    return(RList)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
