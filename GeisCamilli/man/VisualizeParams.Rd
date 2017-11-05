\name{VisualizeParams}
\alias{VisualizeParams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualize Generated Parameter Distributions
}
\description{
This function will show densities of the important parameters creating a generated response pattern.
}
\usage{
VisualizeParams(struct=NA,samples=1000) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{struct}{
  list; These are the settings for the generation of a response pattern.
  }
  \item{samples}{
  numeric; This is the number of samples to use to visualize the distribution sampled.
  }
}
\details{
This is for pedagogical purposes. The generalized ICC for the model fit in this code has the form: 
\newline 
P(Y|Theta) = c + (1-c) * Probit(A*Theta - (b+tau)).
}
\value{
A plot of the distributions of A, b, tau, c, and Theta.
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
VisualizeParams()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
