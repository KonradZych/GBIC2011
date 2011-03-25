\name{genotypesToCross}
\alias{genotypesToCross}

\title{Cross object from genotypic data.}

\description{
  Produces from genotypic matrix file containing object of type cross, reads it into R a returns.
}

\usage{
genotypesToCross(genotypeMatrix, expressionMatrix, doClustering=FALSE, groups=10, iterations = 100, outputFile="mycross.csv", verbose=FALSE, debugMode=0)
}

\arguments{
 \item{genotypeMatrix}{ Matrix with genotype values with: columns, individuals and on the rows, markers }
 \item{expressionMatrix}{ Matrix with expression values with: columns, individuals and on the rows, markers }
 \item{doClustering}{ Cluster the chromosomes into groups }
 \item{groups}{ Number of expected chromosomes }
 \item{iterations}{ how many times clustering should be performed }
 \item{outputFile}{ Filename to save the created cross object to }
 \item{verbose}{ Be verbose}
 \item{debugMode}{ 1: Print our checks, 2: print additional time information }
}

\value{
  A matrix with genotypes
}

\details{
  TODO
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
	Under tender patronage of: Danny Arends \email{Danny.Arends@gmail.com}
}

\examples{
	library(qtl)
  	data(genotypes)
	data(expression_ratio)
	genotypesToCross(genotypes,expression_ratio)
}

\seealso{
  TODO
}

\keyword{manip}
