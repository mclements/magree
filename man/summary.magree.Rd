\name{summary.magree}
\alias{summary.magree}
\alias{summary.oconnell}
\alias{summary.schouten}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{summary method for magree objects
%%  ~~function to do ... ~~
}
\description{summary method for magree objects
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{summary}{magree}(object, ...)
\method{summary}{oconnell}(object, ci.transform = c("logit", "identity"), ci.p = 0.95, ...)
\method{summary}{schouten}(object, ci.transform = c("logit", "identity"), ci.p = 0.95, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{obj} here~~
}
  \item{ci.transform}{transformation used to calculate the confidence
    intervals. Either "logit" for a logit transform or "identity" for no transform.
%%     ~~Describe \code{ci.transform} here~~
}
  \item{ci.p}{p value for the confidence interval.
%%     ~~Describe \code{ci.p} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
%% \details{
%% %%  ~~ If necessary, more details than the description above ~~
%% }
%% \value{
%% %%  ~Describe the value returned
%% %%  If it is a LIST, use
%% %%  \item{comp1 }{Description of 'comp1'}
%% %%  \item{comp2 }{Description of 'comp2'}
%% %% ...
%% }
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
%% \author{
%% %%  ~~who you are~~
%% }
%% \note{
%% %%  ~~further notes~~
%% }

%% %% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% %% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
summary(magree(landis))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ design }% use one of  RShowDoc("KEYWORDS")
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
