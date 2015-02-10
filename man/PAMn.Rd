\name{PAMn}
\alias{PAMn}
\alias{Dayhoffn}

\title{
  PAM and Dayhoff matrices calculation
}
\description{
  These functions calculate the mutation matrix or Dayhoff matrix
  from the mutation matrix at PAM 1 and base backrgound frequency.
}
\usage{
  PAMn(PAM1, n)
  Dayhoffn(PAM1, BF, n)
}

\arguments{
  \item{PAM1}{
  A matrix of \link{numeric}: the mutation probability 
  from one AA to another AA at PAM distance 1.
  The order of AA in the matrix should be identical to \link{AACharacterSet}.
  }
  \item{n}{
    A \link{numeric}: the PAM distance.
  }
  \item{BF}{
    A \link{numeric} vector: the backrgound frequency of AAs.
    The order of AA in the vector should also be identical to 
    \link{AACharacterSet}.
  }
}
\details{
  Calculate the n-PAM matrices from PAM1 mutation matrix and n.
  To compute n-PAM matrices, we multiply the PAM1 matrix through itself
  N times, which is most efficiently achieved through
  n additions in log space.
  
  Computing Dayhoff matrices from PAM mutation matrices and AA frequency.
  Dayhoff matrices are the ratios
  P("alignment i and j arose through evolution") /  P("alignment i and j arose by chance")
}
\value{
  A \link{numeric} matrix is returned.
}
\references{
Dayhoff, M.O., and Schwartz, R.M. (1978). A model of evolutionary change in proteins. In In Atlas of Protein Sequence and Structure,.

Gonnet, G.H., and Scholl, R. (2009). Scientific Computation (Cambridge, UK; New York: Cambridge University Press).

}

\author{
  Ge Tan
}

\examples{
  data(GONNET)
  data(GONNETBF)
  ## PAM 250 mutation matrix
  PAM250 <- PAMn(GONNET, 250)
  
  ## Dayhoff 250 matrix
  Dayhoff250 <- Dayhoffn(GONNET, GONNETBF, 250)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
