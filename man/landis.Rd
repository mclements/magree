\name{landis}
\alias{landis}
\docType{data}
\title{
Landis and Koch dataset.
}
\description{Canonical dataset for agreement for multiple observers described in
  Landis and Koch (Biometrics 1977; 33: 363-374).
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("landis")}
\format{
  The format is:
 int [1:118, 1:7] 4 1 3 4 3 2 1 3 2 1 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:118] "1" "2" "3" "4" ...
  ..$ : chr [1:7] "A" "B" "C" "D" ...
}
%% \details{, described in )
%% %%  ~~ If necessary, more details than the __description__ above ~~
%% }
\source{
Landis and Koch (Biometrics 1977; 33: 363-374)
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
%% \references{
%% %%  ~~ possibly secondary sources and usages ~~
%% }
\examples{
data(landis)
## maybe str(landis) ; plot(landis) ...
}
\keyword{datasets}
