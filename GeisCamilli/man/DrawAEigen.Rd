\name{DrawAEigen}
\alias{DrawAEigen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Eigenanalysis to Sample Slopes (Polytomous)
}
\description{
Eigenanalysis for polytomous slopes with no check on sign.
}
\usage{
DrawAEigen(covZ-diag(J),Q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covZ}{Augmented Data covariance matrix minus diag(J)}
  \item{Q}{Dimension of the latent factor}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a list of 
Atemp = Eigenvectors 
Avec = Eigenvalues
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
