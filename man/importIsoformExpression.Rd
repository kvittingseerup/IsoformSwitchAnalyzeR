\name{importIsoformExpression}
\alias{importIsoformExpression}

\title{
Import expression data from Kallisto, Salmon or RSEM into R.
}

\description{
A general-purpose import function which imports isoform expression data from Kallisto, Salmon or RSEM into R. This is a wrapper for the tximport package with some extra functionalities and is meant to be used to import the data and afterwards a switchAnalyzeRlist can be created with \code{importRdata}. It is highly reccomended that both the imported TxPM and counts values are used both in the creation of the switchAnalyzeRlist with \code{importRdata} (through the "isoformCountMatrix" and "isoformRepExpression" arguments). Most notably this import function also enables inter-library normalization (via edgeR) of the abundance estimates. Note that the pattern argument allows import of only a subset of files.
}

\usage{
importIsoformExpression(
    parentDir,
    calculateCountsFromAbundance=TRUE,
    interLibNormTxPM=TRUE,
    normalizationMethod='TMM',
    pattern='',
    invertPattern=FALSE,
    ignore.case=FALSE,
    showProgress = TRUE,
    quiet = FALSE
)
}

\arguments{
\item{parentDir}{ Parrent directory where each quantified sample is in a sub-directory.}


\item{calculateCountsFromAbundance}{A logic indicating whether to generate estimated counts using the estimated abundances. Recomended as it will incooperate the bias correction algorithms into the analysis. Default is TRUE.}

\item{interLibNormTxPM}{A logic indicating whether to apply an inter-library normalization (via edgeR) to the imported abundances. Recomended as it allow comparison of abundances between samples. Default is TRUE.}

\item{normalizationMethod}{A string indicating the method used for the inter-library normalization. Must be one of "TMM", "RLE", "upperquartile". See \code{?edgeR::calcNormFactors} for more details. Default is "TMM".}

\item{pattern}{ Character string containing a regular expression for which files to import (applied to full path). Defailt is "" corresponding to all. See base::grepl for more details.}

\item{invertPattern}{ Logical. If TRUE return indices or values for elements that do not match..}

\item{ignore.case}{ if FALSE, the pattern matching is case sensitive and if TRUE, case is ignored during matching.}

\item{showProgress}{ A logic indicating whether to make a progress bar (if TRUE) or not (if FALSE). Default is FALSE.}

\item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}
}

\details{
This function requires all data that should be imported is in a directory (as indicated by \code{parentDir}) where each quantified sample is in a seperate sub-directory.

For Kallisto the bias estimation is enabled by adding '--bias' to the function call. For Salmon the bias estimation is enabled by adding '--seqBias' and '--gcBias' to the function call. For RSEM the bias estimation is enabled by adding '--estimate-rspd' to the function call.
}

\value{
A \code{list} containing an abundance matrix, a count matrix and a matrix with the effective lengths for each isoform quantified (rows) in each sample (col) where the first column contains the isoform_ids. The options used for import are stored under the "importOptions" entry). The abundance estimates are in the unit of Transcripts Per Million (TPM) and measrung the relative abundance of a specic transcript.

Transcripts Per Million values are abbreviated to TPM by RSEM, Kallisto and Salmon but will here refered to as TxPM to avoid confusion with the commenly used Tags Per Million (which have been around for way longer). TxPM is an equivilent to RPKM/FPKM except it has been adjusted for as all the biases being modeled by the tools used for the quantifictation including the fragment length distribution and sequence-specific bias as well as GC-fragment bias (this is specific to each tool and how it was run so you need to look up the specific tool). The TxPM is optimal for expression comparison of abundances since most biases will be taking into account.
}

\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
Soneson et al. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research 4, 1521 (2015).
Robinson et al. A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biology (2010)
}
\author{
Kristoffer Vitting-Seerup
}

\seealso{
\code{\link{importRdata}}\cr
\code{\link{createSwitchAnalyzeRlist}}\cr
\code{\link{preFilter}}
}

