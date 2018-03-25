\name{DrawALowerDiag}
\alias{DrawALowerDiag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Lower Triangular Slope Calculation (Polytomous)
}
\description{
Calculate Slopes using a lower triangular matrix constrained to prevent rotation.
}
\usage{
DrawALowerDiag(covT,covTZ,Q,m,n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covT}{Latent factor covariance}
  \item{covTZ}{Covariance of latent factors crossed on augmented Z's}
  \item{Q}{Dimension of the latent factor}
  \item{m}{Number of Items}
  \item{n}{Number of examinees}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a list of 
Atemp = Slope Vectors 
Avec = NA (because it's not an eigenanalysis)
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
