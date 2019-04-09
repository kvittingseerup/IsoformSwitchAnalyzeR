\name{isoformToIsoformFraction}
\alias{isoformToIsoformFraction}
\title{
Calculate isoform fraction from isoform abundance matrix
}
\description{
General purpose function to calculate isoform fraction (IF) matrix from isoform abundace (and potentially gene abundance) matrix.
}
\usage{
isoformToIsoformFraction(
    isoformRepExpression,
    geneRepExpression=NULL,
    isoformGeneAnnotation=NULL,
    quiet = FALSE
)
}

\arguments{
\item{isoformRepExpression}{
A replicate isoform abundace matrix (not log-transformed). The isoform:gene relationship can be provided by either:
\itemize{
    \item{Having \code{isoformRepExpression} contain two additional colmns 'isoform_id' and 'gene_id' indicating which isoforms are a part of which gene}
    \item{Using the {isoformGeneAnnotation} argument}.
}
Importantly \code{isoformRepExpression} must contain isoform ids either as seperate colum called 'isoform_id' or as row.names. The function will figure it out by itself in what combination the annotation is supplied.
}

\item{geneRepExpression}{
Optional. A gene replciate abundance matrix. Must contain gene ids either as seperate colum called 'gene_id' or as row.names.
}

\item{isoformGeneAnnotation}{
A data.frame or GRange with two (meta) collumns: 'isoform_id' and 'gene_id' indicating the relationship between isoforms and parent gene.
}

\item{quiet}{ A logic indicating whether to avoid printing progress messages. Default is FALSE}
}

\details{
This function calculates isoform fractions from isoform abundances. If \code{geneRepExpression} is not supplied the function automatically calculate it by itself.

Note that:
1) isoform:gene relationship can be supplied as two columns either in the \code{isoformRepExpression} or as a seperate data.frame to \code{isoformGeneAnnotation}.
2) The ids in \code{isoformRepExpression} and \code{geneRepExpression} can be supplied either as row.names or as seperate colums respectively called 'isoform_id' and 'gene_id'.
}

\value{
A replicate isoform fraction matrix with layout similar to \code{isoformRepExpression}
}

\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
}
\author{
Kristoffer Vitting-Seerup
}

\examples{
### Please note
# 1) The way of importing files in the following example with
#       "system.file('pathToFile', package="IsoformSwitchAnalyzeR") is
#       specialiced to access the sample data in the IsoformSwitchAnalyzeR package
#       and not somhting you need to do - just supply the string e.g.
#       "myAnnotation/isoformsQuantified.gtf" to the functions
# 2) importRdata directly supports import of a GTF file - just supply the
#       path (e.g. "myAnnotation/isoformsQuantified.gtf") to the isoformExonAnnoation argument

### Import quantifications
salmonQuant <- importIsoformExpression(system.file("extdata/", package="IsoformSwitchAnalyzeR"))


### Extract gene info
localAnnotaion <- rtracklayer::import(system.file("extdata/example.gtf.gz", package="IsoformSwitchAnalyzeR"))[,c('transcript_id','gene_id')]
colnames(localAnnotaion@elementMetadata)[1] <- 'isoform_id'

### Calculate isoform fractions
repIF <- isoformToIsoformFraction(
    isoformRepExpression  = salmonQuant$abundance,
    isoformGeneAnnotation = localAnnotaion
)
}
