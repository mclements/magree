\name{schouten}
\alias{schouten}
\title{
Schouten estimators for multirater agreement.
}
\description{
  Use the Schouten estimator of agreement for nominal or ordinal
  data. This includes a range of statistics on agreement. 
}
\usage{
schouten(X, weights=c("unweighted","linear","quadratic","user"), w=NULL)
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
      agreement weights \eqn{w_{i,j}=1-|i-j|/(c-1)}.}
    \item{"quadratic"}{For ordinal categories where disagreement is proportional to
      the square of the distance between the categories. This is analogous
      to the agreement weights \eqn{w_{i,j}=1-(i-j)^2/(c-1)^2}.}
    \item{"user"}{An indicator for a user-defined weight matrix. The
      \code{weights} argument will be defined as "user" if the \code{w}
      argument is specified.}
  }
}
\item{w}{A user-defined weights matrix. This argument takes precedence
  over \code{weights} if it is specified and the \code{weight} argument
  will be defined as "user".
}
}
\details{
  Fortran code was written based on Schouten (1982).
  
  The output object is very similar to the Fortan code. Not
  all of the variance terms are currently used in the \code{print} and \code{summary} methods.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{N}{Number of subjects (items)}
\item{M}{Number of observers (raters)}
\item{L}{Number of categories}
\item{data}{Re-formatted X}
\item{w}{Weight matrix}
\item{kab}{Kappas between each pair of raters}
\item{ka}{Average kappas for each rater}
\item{kappa}{Average kappa}
\item{pab,pa,p,ma,qab,qa,q,oab,eab,oa,ea,o,e,wa,wab}{Working fields}
\item{varkab}{Variances for kab}
\item{varka}{Variances for ka}
\item{vark}{Variance for the kappa}
\item{covkka}{Covariance term between the overall average kappa and the
  average kappas for each rater}
\item{chi}{Chi-squared statistics comparing the overall average kappa
  and the average kappa for each rater (df=1 under the null hypothesis)}
\item{pchi}{P-values that the overall average kappa
  equals the average kappa for each rater}
\item{var0kab}{Variance for kab under the null hypothesis}
\item{var0ka}{Variance for ka under the null hypothesis}
\item{var0k}{Variance for the overall average kappa under the null hypothesis}
\item{p0}{P-value for kappa=0}
\item{p0a}{P-values that the average kappa for a rater equals zero
  (i.e. ka=0)}
\item{weights}{As input}
\item{X}{As input}
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

## Weights matrix used by Schouten (1982)
w <- outer(1:5,1:5,function(x,y) ((x<=2 & y<=2) | (x>=3 & y>=3))+0)
fit <- schouten(landis,w=w) # user-defined weights

summary(fit) # Schouten (1982), Tables 2 and 5
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
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
