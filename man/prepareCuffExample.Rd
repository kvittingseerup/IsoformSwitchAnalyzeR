\name{prepareCuffExample}
\alias{prepareCuffExample}
\title{Prepare the Cufflinks example data}
\description{
Prepare the Cufflinks example data set.
}
\usage{
prepareCuffExample()
}

\details{
This helper function prepares the Cufflinks example dataset, including the example GTF-file.
}

\value{
A CuffSet object.
}

\author{Kristoffer Vitting-Seerup, Johannes Waage}

\references{
Vitting-Seerup K , Porse BT, Sandelin A, Waage J. (2014) spliceR: an R package for classification of alternative splicing and prediction of coding potential from RNA-seq data. BMC Bioinformatics 15:81.
}


\examples{
#Load cufflinks example data
cuffDB <- prepareCuffExample()
}