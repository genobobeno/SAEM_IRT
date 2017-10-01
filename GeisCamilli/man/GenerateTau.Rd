\name{GenerateTau}
\alias{GenerateTau}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
Generate the offsets from the intercept for polythomous responses.
}
\usage{
GenerateTau(j,ncat,taudist,tauparams)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{j}{
numeric; Number of items
}
  \item{ncat}{
numeric; Number of categories 
}
  \item{taudist}{
character; Distribution of the tau intercepts
}
  \item{tauparams}{
numeric vector; Parameters of the distribution of the tau intercepts
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
vector of numeric values}
\references{
Geis, Camilli
}
\author{
Eugene Geis
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{GenerateRP}}, \code{\link{GenerateTheta}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
