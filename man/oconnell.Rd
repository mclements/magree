\name{oconnell}
\alias{oconnell}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Implement the O'Connell-Dobson estimator for multirater agreement.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
oconnell(X, score = NULL, i = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A matrix or data-frame with observations/subjects as rows and raters as columns.
}
  \item{score}{
    The scores that are to be assigned to the categories. Currently,
    this defaults to a sorted list of the unique values.
}
\item{i}{
  \enumerate{
    \item{For nominal categories - only perfect agreement is counted.}
    \item{For ordinal categories where disagreement is proportional to
      the distance between the categories. This is analogous to the
      agreement weights \eqn{w_{i,j}=1-|i-j|/(c-1)}.}
    \item{For ordinal categoreis where disagreement is proportional to
      the square of the distance between the categories. This is analogous
      to the agreement weights \eqn{w_{i,j}=1-(i-j)^2/(c-1)^2}.}
  }
  }
}
\details{
The Fortran code from Professor Dianne O'Connell was adapted for R. 

Currently, the output object is very similar to the Fortan code.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{nrater}{Number of observers}
\item{nscore}{Number of categories}
\item{nsubj}{Number of subjects (items)}
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
\item{vars1}{Unconstrained variance of S when observers are homogeneous}
\item{v0sav1}{Null variance of Sav when observers are distinct}
\item{v0sav2}{Null variance of Sav when observers are homogeneous}
\item{vsav1}{Unconstrained variance of Sav when observers are distinct}
\item{vsav2}{Unconstrained variance of Sav when observers are homogeneous}
%%\item{resp[i,j]}{Response for observer i on subject j; transpose of X (BEWARE)}
\item{X[j,i]}{Response for observer i on subject j}
\item{score(i)}{Score associated with i'th category}
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

fit <- oconnell(landis, i=1)
fit3 <- oconnell(landis, i=3)
## Table 1 (O'Connell and Dobson, 1984)
summary(fit)

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
## table(fit$d)
## tapply(names(fit$d),fit$d,I)
data.frame(SlideType=1:18,
           S1=sapply(slideTypeGroups, function(ids) mean(fit$s1[as.character(ids)])),
           S2=sapply(slideTypeGroups, function(ids) mean(fit$s2[as.character(ids)])))

## Table 5, O'Connell and Dobson (1984)
oconnell(landis==1)
oconnell(landis==2)
oconnell(landis==3)
oconnell(landis==4)
oconnell(landis==5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
