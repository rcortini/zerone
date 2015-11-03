\name{zerone}
\alias{zerone}
\title{Zerone}
\description{
   Discretize multiple ChIP-seq profiles while performing a quality control.
}
\usage{
zerone(y)
}
\arguments{
    \item{y}{data.frame with at least three columns containing the data of
    every window in your profiles. (See Details)}
}
\details{
    The first column of \var{y} must contain the chromosome name of each
    window, the second column must contain the summed read count of all your
    control profiles, and the third and subsequent columns the read count of
    each of your target profiles.
}
\value{
    A numeric vector containing the Viterbi path: the most likely state of
    each window, where \bold{State 2 represents an enrichment}.
}
\author{
    Pol Cuscó and Guillaume Filion. date: October 29, 2015
}
\references{
    Cuscó, P. and Filion, G. (2015) Zerone: a ChIP-seq discretizer with
    built-in quality control. Journal. \bold{vol},page--page
}
\examples{
## Run this from the zerone directory:

library(zerone)
path <- zerone(ZeroneExampleData)

## The window 534 is enriched (has state 2 in the Viterbi path).

path[530:540]
# [1] 1 1 1 1 2 0 0 1 1 1 1
}
