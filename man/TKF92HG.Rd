\name{TKF92HG}
\alias{TKF92HG}
\alias{TKF92HGPair}

\title{
  Evolutionary distance estimation with TKF92 model considering the regional heterogeneity of substitution rates
}
\description{
  This function implements the TKF92 model to estimate the pairwise distance
  from protein sequences.
  An additional simple model of regional heterogeneity of substitution rates is used.
}
\usage{
TKF92HG(fasta, mu=NULL, r=NULL, Ps=NULL, Kf=NULL, expectedLength=362, 
        substModel, substModelBF)
TKF92HGPair(seq1, seq2, mu=NULL, r=NULL, Ps=NULL, Kf=NULL, distance=NULL,
            expectedLength=362, substModel, substModelBF)
}

\arguments{
  \item{fasta}{
    A named list of sequences in vector of characters format.
  \code{read.fasta} from package \code{seqinr} outputs this format 
  when reading from a fasta file.
  }
  \item{mu}{
    A numeric value between 0 and 1 or NULL.
    It is the death rate per normal link in TKF92 model.
    When it is NULL, a joint estimation of \code{mu}, \code{r}, \code{Ps}, \code{Kf} 
    and \code{distance} will be done.
    When it is given, only the distance will be estimated.
  }
  \item{r}{
    A numeric value between 0 and 1 or NULL.
    It is the success probability of the geometric distribution 
    for modeling the fragment length in TKF92 model.
    When it is NULL, a joint estimation of \code{mu}, \code{r}, \code{Ps}, \code{Kf} 
    and \code{distance} will be done.
    When it is given, only the distance will be estimated.
  }
  \item{Ps}{
    A numeric value between 0 and 1 or NULL.
    It is the equilibrium frequency of slow fragments.
    Hence the equilibrium frequency of fast fragements is Pf= 1- Ps.
    When it is NULL, a joint estimation of \code{mu}, \code{r}, \code{Ps}, \code{Kf} 
    and \code{distance} will be done.
    When it is given, only the distance will be estimated.
  }
  \item{Kf}{
    A numeric value larger than 1 or NULL.
    It is the ratio of substitutions rates between fast fragments and slow fragments.
    When it is NULL, a joint estimation of \code{mu}, \code{r}, \code{Ps}, \code{Kf} 
    and \code{distance} will be done.
    When it is given, only the distance will be estimated.
  }
  \item{distance}{
    A numeric value: the PAM distance between two protein sequences.
    When it is given, TKF92HGPair only calculates the negative log-likelihood.
  }
  \item{expectedLength}{
    A numeric object: the expected length of input protein sequences.
    By default, the average sequence length, 362, from OMA browser is used.
  }
  \item{substModel}{
    A numeric matrix: the mutation probability from one AA
    to another AA at PAM distance 1. 
    The order of AA in the matrix should be identical to \link{AACharacterSet}.
  }
  \item{substModelBF}{
    A vector of numeric: the backrgound frequency of AAs.
    The order of AA in the vector should also be identical to 
    \link{AACharacterSet}.
  }
  \item{seq1, seq2}{
    A vector of character: the sequences of two proteins to compare.
  }
}

\details{
  Currently this implementation only supports the normal 20 AAs.
  Missing or Ambiguous characters are not supported.
  
  This is a very simple model of substitution rate heterogeneity.
  This model assumes that there are only two varieties of fragments: 
  one with relatively fast substitution rates and 
  the other with slow substitution rates.
  This model also assumes the fragment size distribution of these two fragments is identical.
}

\value{
  A list of matrices are returned: the matrix of estimated distances, 
  the matrix of estimated distance variances, 
  the matrix of negative log-likelihood between the sequences.
}

\references{
Thorne, J.L., Kishino, H., and Felsenstein, J. (1992). Inching toward reality: an improved likelihood model of sequence evolution. J. Mol. Evol. 34, 3-16.

Gonnet, G.H., Cohen, M.A., and Benner, S.A. (1992). Exhaustive matching of the entire protein sequence database. Science 256, 1443-1445.
}
\author{
  Ge Tan
}


\seealso{
  \code{\link{AACharacterSet}},
  \code{\link{GONNET}}, \code{\link{GONNETBF}}
}
\examples{
  \donttest{
    ## This example is not tested due to running time > 5s
  data(GONNET)
  data(GONNETBF)
  library(seqinr)
  fasta <- read.fasta(file.path(system.file("extdata", package="TKF"),
                      "pair1.fasta"),
                      seqtype="AA", set.attributes=FALSE)
  ## 1D estimation: only distance
  TKF92HG(fasta, mu=5.920655e-04, r=0.8, Ps=1, Kf=1.2,
          substModel=GONNET, substModelBF=GONNETBF)
  
  ## 2D estimation: joint estimation of distance, mu and r
  TKF92HG(fasta, substModel=GONNET, substModelBF=GONNETBF)
  
  ## only apply to a pair of sequences
  seq1 <- fasta[[1]]
  seq2 <- fasta[[2]]
  TKF92HGPair(seq1, seq2, mu=5.920655e-04, r=0.8, Ps=1, Kf=1.2,
            substModel=GONNET, substModelBF=GONNETBF)
  }
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
