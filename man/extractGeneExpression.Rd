\name{extractGeneExpression}
\alias{extractGeneExpression}
\title{
Extract raw gene counts or abundances from a switchAnalyzeRlist object.
}
\description{
Extract replicate gene raw unnormalised counts or expression from a switchAnalyzeRlist object using all the annotation fixes employed in creating the switchAnalyzeRlist.
}
\usage{
extractGeneExpression(
    switchAnalyzeRlist,
    extractCounts = TRUE,
    addGeneNames = TRUE,
    addIdsAsColumns = TRUE
)
}

\arguments{
\item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object.
}

\item{extractCounts}{
A logic to indicate whether to extract raw unnormalised counts (if TRUE, default) or expression estimates (if FALSE).
}

\item{addGeneNames}{
A logic to indicate whether to add gene_names to the expression matrix (if TRUE, default) or not (if FALSE).
}

\item{addIdsAsColumns}{
A logic to indicate whether to add the gene identifiers to the data.frame as collumns (if TRUE, default) or rownames (if FALSE).
}

}

\details{
The count matrix obtained if \code{extractCounts=TRUE} is the same as would be obtained by running tximport with \code{countsFromAbundance="scaledTPM" which are suitable both for analysis of differential expression and usage.}
}

\value{
A data.frame with the replicate count/abundance estimates as well as gene_id (and gene_name if \code{extractCounts=TRUE})
}

\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
}
\author{
Kristoffer Vitting-Seerup
}

\examples{
data("exampleSwitchList")

### Raw count matrix
geneCountMatrix <- extractGeneExpression(
    exampleSwitchList,
    extractCounts = TRUE
)

### Raw count matrix - with ids as rownames instead of columns
geneCountMatrix <- extractGeneExpression(
    exampleSwitchList,
    extractCounts = TRUE,
    addIdsAsColumns = FALSE
)

### Abundance matrix
geneExpresionMatrix <- extractGeneExpression(
    exampleSwitchList,
    extractCounts = FALSE
)
}
