\name{print.oconnell}
\alias{print.oconnell}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print out a short description of the \code{oconnell} object.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.oconnell(X, score = NULL, i = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{score}{
%%     ~~Describe \code{score} here~~
}
  \item{i}{
%%     ~~Describe \code{i} here~~
}
}
\details{
The print report is currently very short.
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
function (X, score = NULL, i = 3) 
{
    nrater <- ncol(X)
    nsubj <- nrow(X)
    X <- matrix(as.integer(as.matrix(X)), nrow(X), ncol(X))
    if (is.null(score)) 
        score <- sort(unique(as.vector(X)))
    nscore <- length(score)
    .Fortran("oconnell", nrater = nrater, nsubj = nsubj, nscore = nscore, 
        i = as.integer(i), r = rep(0L, nrater), resp = t(X), 
        s1 = rep(0, nsubj), s2 = rep(0, nsubj), p1 = matrix(0, 
            nrater, nscore), p2 = rep(0, nscore), w1 = matrix(0, 
            nrater, nscore), w2 = matrix(0, nrater, nscore), 
        score = as.double(score), delta = matrix(0, nrater, nrater), 
        d = rep(0, nsubj), expd1 = 0, expd2 = 0, dbar = 0, vars1 = 0, 
        var0s1 = 0, vsav1 = 0, v0sav1 = 0, vars2 = 0, var0s2 = 0, 
        vsav2 = 0, v0sav2 = 0, sav1 = 0, sav2 = 0)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
