\name{QueryUser}
\alias{QueryUser}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Query User for Input 
}
\description{
Used to validate user input or to request more input
}
\usage{
QueryUser(message,choices=NA,type="character",defaultchoice=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{message}{
character; Tell the user what you want
}
  \item{choices}{
vector; populate choices or leave empty (NA) for open-ended
}
  \item{type}{
character; "numeric" or "character"
}
  \item{defaultchoice}{
numeric or character; dependent on type
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
if type=="character" and choices is a vector, user's input number will return that indexed element in the choices vector
if type=="character" and choices=NA, it will return user's input string
if type=="numeric" and choices=NA, it will return user's input number
if type=="numeric" and choices is a vector, user's input number will return that indexed element in the choices vector
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
QueryUser("Enter a number:", choices=c(4,3,2,1),type="numeric",defaultchoice=4)
QueryUser("Enter a number:", choices=NA,type="numeric",defaultchoice=4)
QueryUser("Enter a number:", choices=c(4,3,2,1),type="character",defaultchoice=4)
QueryUser("Enter a number:", choices=NA,type="character",defaultchoice=4)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
