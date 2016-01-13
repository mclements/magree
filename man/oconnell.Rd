\name{oconnell}
\alias{oconnell}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
O'Connell-Dobson estimators for multiobserver agreement.
}
\description{
  Use the O'Connell-Dobson estimator of agreement for nominal or ordinal
  data. This includes a range of statistics on agreement for assuming
  either distinct or homogeneous items. 
}
\usage{
oconnell(X, weights=c("unweighted","linear","quadratic"), i=NULL, score = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A matrix or data-frame with observations/subjects as rows and observers as columns.
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
\item{i}{
  \enumerate{
    \item{For nominal categories - only perfect agreement is counted.}
    \item{For ordinal categories where disagreement is proportional to
      the distance between the categories. This is analogous to the
      agreement weights \eqn{w_{i,j}=1-|score[i]-score[j]|/(max(score)-min(score))}.}
    \item{For ordinal categories where disagreement is proportional to
      the square of the distance between the categories. This is analogous
      to the agreement weights \eqn{w_{i,j}=1-(score[i]-score[j])^2/(max(score)-min(score))^2}.}
  }
  This argument takes precedence over \code{weights} if it is specified.
}
  \item{score}{
    The scores that are to be assigned to the categories. Currently,
    this defaults to \code{1:L}, where code{L} is the number of categories.
}
}
\details{
  The Fortran code from Professor Dianne O'Connell was adapted for R. 
  
  The output object is very similar to the Fortan code. Not
  all of the variance terms are currently used in the \code{print},
  \code{summary} and \code{plot} methods.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{X}{As input}
\item{i}{As input}
\item{nrater}{Number of observers}
\item{nscore}{Number of categories}
\item{nsubj}{Number of subjects}
\item{p1[j,k]}{Probability of observer j giving score k when observers
  are distinct}
\item{p2[k]}{Probability of score k when observers are homogeneous}
\item{w1[j,k]}{Weighted average of d[] for observer j, score k}
\item{w2[k]}{Weighted average of d[] for score k when observers are homogeneous}
\item{d[j]}{Amount of disagreement for subject j}
\item{s1[j]}{Chance-corrected agreement statistic for subject j when
  observers are distinct}
\item{s2[j]}{Chance-corrected agreement statistic for subject j when
  observers are homogeneous; s[j]=1-d[j]/expdel.}
\item{delta[j,k]}{j<k: amount of disagreement expected by change for
  observers j and k; j>k amount of disagreement expected by chance for
  observers j and k when observers are homogeneous}
\item{expd1}{Amount of disagreement expected by chance in null case when
observers are distinct}
\item{expd2}{Amount of disagreement expected by chance when observers
  are homogeneous}
\item{dbar}{Average value of d[] over all subjects}
\item{sav1}{Chance-corrected agreement statistic over all subjects when
  observers are distinct}
\item{sav2}{Chance-corrected agreement statistic over all subjects when
  observers are homogeneous}
\item{var0s1}{Null variance of S when observers are distinct}
\item{var0s2}{Null variance of S when observers are homogeneous}
\item{vars1}{Unconstrained variance of S when observers are distinct}
\item{vars2}{Unconstrained variance of S when observers are homogeneous}
\item{v0sav1}{Null variance of Sav when observers are distinct}
\item{v0sav2}{Null variance of Sav when observers are homogeneous}
\item{vsav1}{Unconstrained variance of Sav when observers are distinct}
\item{vsav2}{Unconstrained variance of Sav when observers are homogeneous}
\item{p0sav1}{Probability of overall agreement due to chance when
  observers are distinct}
\item{p0sav2}{Probability of overall agreement due to chance when
  observers are homogeneous}
\item{resp[i,j]}{Response for observer i on subject j; transpose of X (BEWARE)}
\item{score(i)}{Score associated with i'th category}
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
\code{\link{magree}}, \code{\link{schouten}}.
}
\examples{

## Table 1 (O'Connell and Dobson, 1984)
summary(fit <- oconnell(landis, weights="unweighted"))
update(fit, weights="linear")
update(fit, weights="quadratic")

## Table 3 (O'Connell and Dobson, 1984)
slideTypeGroups <-
    list(c(2,3,5,26,31,34,42,58,59,67,70,81,103,120),
         c(7,10:13,17,23,30,41,51,55,56,60,65,71,73,76,86,87,105,111,116,119,124),
         c(4,6,24,25,27,29,39,48,68,77,79,94,101,102,117),
         c(9,32,36,44,52,62,84,95),
         c(35,53,69,72),
         c(8,15,18,19,47,64,82,93,98,99,107,110,112,115,121),
         c(1,16,22,49,63,66,78,90,100,113),
         c(28,37,40,61,108,114,118),
         106,
         43,
         83,
         c(54,57,88,91,126),
         c(74,104),
         38,
         46,
         c(89,122),
         c(80,92,96,123),
         85)
data.frame(SlideType=1:18,
           S1=sapply(slideTypeGroups,
                     function(ids) mean(fit$s1[as.character(ids)])),
           S2=sapply(slideTypeGroups,
                     function(ids) mean(fit$s2[as.character(ids)])))

## Table 5, O'Connell and Dobson (1984)
oconnell(landis==1)
oconnell(landis==2)
oconnell(landis==3)
oconnell(landis==4)
oconnell(landis==5)

## Plot of the marginal distributions
plot(fit)
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
