\name{WrapT}
\alias{WrapT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Sampler for Latent Factor (Polytomous)
}
\description{
Function samples one examinee from multivariate distribution for Latent Factor conditional on current parameters.
}
\usage{
WrapT(i,A,Z,BTB_INV,b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{i}{
Examinee's row index
}
  \item{A}{Current Latent Factor Loadings
}
  \item{Z}{Augmented Data Matrix}
  \item{BTB_INV}{Normalization Matrix}
  \item{b}{Current Intercepts}
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
