\name{AAToInt}
\alias{AAToInt}
\alias{AACharacterSet}
\alias{AAGapCharacterSet}
\alias{AmbiguousAACharacterSet}
\alias{AmbiguousAAGapCharacterSet}
\alias{DNACharacterSet}
\alias{DNAGapCharacterSet}
\alias{RNACharacterSet}
\alias{RNAGapCharacterSet}

\title{
  AA, DNA, RNA character set
}
\description{
  Some AA, DNA, RNA character set defined in this package and 
  functions to convert them into integers.
}
\usage{
AAToInt(AA)
}

\arguments{
  \item{AA}{
    A vector of \link{character}.
  }
}
\details{
  Each AA is converted to the position of that AA in \code{AACharacterSet}.
}

\value{
  A \code{integer} vector.
}

\author{
  Ge Tan
}

\examples{
  library(seqinr)
  fasta <- read.fasta(file.path(system.file("extdata", package="TKF"),
                      "pair1.fasta"),
                      seqtype="AA", set.attributes=FALSE)
  AAToInt(fasta[[1]])
  
  AACharacterSet  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
