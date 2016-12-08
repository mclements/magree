\name{plot.magree}
\alias{plot.magree}
\alias{plot.oconnell}
\alias{plot.schouten}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot methods for magree, oconnell and schouten objects
}
\description{
plot methods for magree, oconnell and schouten objects
}
\usage{
\method{plot}{magree}(x, type = c("p1", "kappa by observer"),
xlab = NULL, ylab = NULL, main = NULL, ...)
\method{plot}{oconnell}(x, type = c("p1"), xlab = NULL, ylab = NULL, main = NULL, ...)
\method{plot}{schouten}(x, type = c("kappa by observer"), xlab = NULL,
ylab = NULL,
main = NULL, xdelta = 0.1, axes = TRUE, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ magree, oconnell or schouten object.
%%     ~~Describe \code{obj} here~~
}
  \item{type}{Type of plot. For \code{"p1"}, plot the probabilities by observer. For \code{"kappa by observer"}, plot the kappas for each observer.
%%     ~~Describe \code{type} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{main}{
%%     ~~Describe \code{main} here~~
}
  \item{xdelta}{For plot.schouten and "kappa by observer", specifies the width of
    the brackets for the confidence intervals.
%%     ~~Describe \code{xdelta} here~~
}
  \item{axes}{
Bool for whether to plot the axes.
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
fit <- schouten(landis)
plot(fit)
fit <- oconnell(landis)
plot(fit,type="p1")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ design }% use one of  RShowDoc("KEYWORDS")
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
