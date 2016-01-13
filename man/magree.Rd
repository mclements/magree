\name{magree}
\alias{magree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
O'Connell-Dobson-Schouten estimators for multiobserver agreement.
}
\description{
  Use the O'Connell-Dobson-Schouten estimators of agreement for nominal or ordinal
  data.  
}
\usage{
magree(X, weights=c("unweighted","linear","quadratic"), score = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A matrix or data-frame with observations/subjects as rows and raters as columns.
}
\item{weights}{
  \describe{
    \item{"unweighted"}{For nominal categories - only perfect agreement is counted.}
    \item{"linear"}{For ordinal categories where disagreement is proportional to
      the distance between the categories. This is analogous to the
      agreement weights \eqn{w_{i,j}=1-|score[i]-score[j]|/(max(score)-min(score))}.}
    \item{"quadratic"}{For ordinal categories where disagreement is proportional to
      the square of the distance between the categories. This is analogous
      to the agreement weights \eqn{w_{i,j}=1-(score[i]-score[j])^2/(max(score)-min(score))^2}.}
  }
}
  \item{score}{
    The scores that are to be assigned to the categories. Currently,
    this defaults to a sorted list of the unique values.
}
}
\details{
  The Fortran code from Professor Dianne O'Connell was adapted for R. 
  
  The output object is very similar to the Fortan code. Not
  all of the variance terms are currently used in the \code{print} and \code{summary} methods.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{oconnell}{object from the \code{oconnell} function}
\item{schouten}{object from the \code{schouten} function}
\item{call}{As per \code{sys.call()}, to allow for using \code{update}}
%% ...
}
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
%% \author{
%% %%  ~~who you are~~
%% }
%% \note{
%% %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{oconnell}}, \code{\link{schouten}}.
}
\examples{

## Table 1 (O'Connell and Dobson, 1984)
summary(fit <- magree(landis, weights="unweighted"))
update(fit, weights="linear")
update(fit, weights="quadratic")

## Table 5, O'Connell and Dobson (1984)
magree(landis==1)
magree(landis==2)
magree(landis==3)
magree(landis==4)
magree(landis==5)

## Plot of the marginal distributions
plot(fit) 

## Plot of the average agreement by observer
plot(fit, type="kappa by observer") 
}
%% ds <- transform(data.frame(id=rownames(landis),
%%                            d1=oconnell(landis,i=1)$s2,
%%                            d2=oconnell(landis,i=2)$s2
%%                            d3=oconnell(landis,i=3)$s2),
%%                 score=d1+d2+d3)
%% ds <- ds[with(ds,order(-d1,d2)),]
%% lapply(unique(ds$score), function(score) ds$id[ds$score==score])
%%
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Agreement }
\keyword{ Kappa }
