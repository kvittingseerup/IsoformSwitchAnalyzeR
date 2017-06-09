\name{extractExpressionMatrix}
\alias{extractExpressionMatrix}
\title{
Extract Gene/Isoform Expression Matrix.
}
\description{
Extract a data.frame with (mean) gene expression, isoform expression or Isoform Fraction values for all conditions (columns) from a switchAnalyzeRlist.
}
\usage{
extractExpressionMatrix(
    switchAnalyzeRlist,
    feature='isoformUsage',
    addInfo=FALSE,
    na.rm=TRUE
)

}

\arguments{
\item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object.
}
  \item{feature}{
The feature of which to extract the expression matrix for. Can be ether 'geneExp' for gene expression levels, 'isoformExp' for isoform expression levels or 'isoformUsage' for IF values. Default is 'isoformUsage'.
}
  \item{addInfo}{
A logic indicating whether annotated non-conditional data (such as gene name, PTC status etc.) should be added to the data.frame. Default is FALSE.
}
  \item{na.rm}{
A logic indicating whether rows with NA expression values should be removed. Default is TRUE.
}

}

\value{
This function returns a data.frame where the first column is the gene/isoform id followed by the mean (if calculated by any of the \code{import*()} functions) expression/usage in all different conditions (one column pr condition) and if addInfo=TRUE then the additional non-conditional dependent data is added as well.
}

\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
}
\author{
Kristoffer Vitting-Seerup
}

\examples{
data("exampleSwitchListAnalyzed")

### Gene expression matrix
geneMatrix <- extractExpressionMatrix(exampleSwitchListAnalyzed, feature = 'geneExp')
head(geneMatrix)

# with additional info
geneMatrix <- extractExpressionMatrix(exampleSwitchListAnalyzed, feature = 'geneExp', addInfo = TRUE)
head(geneMatrix)

### Isoform Fraction value expression matrix
ifMatrix <- extractExpressionMatrix(exampleSwitchListAnalyzed, feature = 'isoformUsage')
head(ifMatrix)

# with additional info
ifMatrix <- extractExpressionMatrix(exampleSwitchListAnalyzed, feature = 'isoformUsage', addInfo = TRUE)
head(ifMatrix)
}
