\name{sechol}
\alias{sechol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sechol(A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
%%     ~~Describe \code{A} here~~
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
function (A) 
{
    n <- nrow(A)
    L <- matrix(rep(0, n * n), ncol = ncol(A))
    macheps <- 2.23e-16
    tau <- macheps^(1/3)
    gamm <- max(A)
    deltaprev <- 0
    Pprod <- diag(n)
    if (n > 2) {
        for (k in 1:(n - 2)) {
            if ((min(diag(A[(k + 1):n, (k + 1):n]) - A[k, (k + 
                1):n]^2/A[k, k]) < tau * gamm) && (min(svd(A[(k + 
                1):n, (k + 1):n])$d)) < 0) {
                dmax <- order(diag(A[k:n, k:n]))[(n - (k - 1))]
                if (A[(k + dmax - 1), (k + dmax - 1)] > A[k, 
                  k]) {
                  print(paste("iteration:", k, "pivot on:", dmax, 
                    "with absolute:", (k + dmax - 1)))
                  P <- diag(n)
                  Ptemp <- P[k, ]
                  P[k, ] <- P[(k + dmax - 1), ]
                  P[(k + dmax - 1), ] = Ptemp
                  A <- P \%*\% A \%*\% P
                  L <- P \%*\% L \%*\% P
                  Pprod <- P \%*\% Pprod
                }
                g <- rep(0, length = (n - (k - 1)))
                for (i in k:n) {
                  if (i == 1) 
                    sum1 <- 0
                  else sum1 <- sum(abs(A[i, k:(i - 1)]))
                  if (i == n) 
                    sum2 <- 0
                  else sum2 <- sum(abs(A[(i + 1):n, i]))
                  g[i - (k - 1)] <- A[i, i] - sum1 - sum2
                }
                gmax <- order(g)[length(g)]
                if (gmax != k) {
                  print(paste("iteration:", k, "gerschgorin pivot on:", 
                    gmax, "with absolute:", (k + gmax - 1)))
                  P <- diag(ncol(A))
                  Ptemp <- P[k, ]
                  P[k, ] <- P[(k + dmax - 1), ]
                  P[(k + dmax - 1), ] = Ptemp
                  A <- P \%*\% A \%*\% P
                  L <- P \%*\% L \%*\% P
                  Pprod <- P \%*\% Pprod
                }
                normj <- sum(abs(A[(k + 1):n, k]))
                delta <- max(0, deltaprev, -A[k, k] + max(normj, 
                  tau * gamm))
                if (delta > 0) {
                  A[k, k] <- A[k, k] + delta
                  deltaprev <- delta
                }
            }
            L[k, k] <- A[k, k] <- sqrt(A[k, k])
            for (i in (k + 1):n) {
                L[i, k] <- A[i, k] <- A[i, k]/L[k, k]
                A[i, (k + 1):i] <- A[i, (k + 1):i] - L[i, k] * 
                  L[(k + 1):i, k]
                if (A[i, i] < 0) 
                  A[i, i] <- 0
            }
        }
    }
    A[(n - 1), n] <- A[n, (n - 1)]
    eigvals <- eigen(A[(n - 1):n, (n - 1):n])$values
    delta <- max(0, deltaprev, -min(eigvals) + tau * max((1/(1 - 
        tau)) * (max(eigvals) - min(eigvals)), gamm))
    if (delta > 0) {
        print(paste("delta:", delta))
        A[(n - 1), (n - 1)] <- A[(n - 1), (n - 1)] + delta
        A[n, n] <- A[n, n] + delta
        deltaprev <- delta
    }
    L[(n - 1), (n - 1)] <- A[(n - 1), (n - 1)] <- sqrt(A[(n - 
        1), (n - 1)])
    L[n, (n - 1)] <- A[n, (n - 1)] <- A[n, (n - 1)]/L[(n - 1), 
        (n - 1)]
    L[n, n] <- A[n, n] <- sqrt(A[n, n] - L[n, (n - 1)]^2)
    return(t(Pprod) \%*\% t(L) \%*\% t(Pprod))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
