\name{zerone}
\alias{zerone}
\title{zerone}
\usage{
zerone(directory)
}
\arguments{
  \item{directory}{directory where the wig files are. ChIP-seq Input profiles must contain the strings \code{input} or \code{control} (case-insensitive) in their file name in order to be recognized as Input.}
}
\value{
  A fitted HMM using the \code{zerone} discretizer. (See \code{zerone})
}
\description{
  Wrapper for the zerone discretizer. Run zerone directly from wig files.
}
\examples{
fit <- wig2zerone("./exampleDirectory")
}
\author{
  Pol Cuscó and Guillaume Filion. date: October 27, 2015.
}
\seealso{
  \code{\link{zerone}}
}
\details{
  The ChIP-seq profiles to discretize must be in the same directory. They consist in one or more Input profiles that will be used as baseline, and one or more Target profiles (these must have been obtained from the same DNA-binding protein and, ideally, they should be experimental replicates). Only \code{fixedStep} formatted wig files are supported.
}
