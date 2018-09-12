\name{removeAnnoationData}
\alias{removeAnnoationData}
\title{
Removes data not needed after switch identification.
}
\description{
This function can removes the biological sequence and/or the replicate quantification data from a switchAnalyzeRlist which will significantly reduce the size of the object. This is intended to be used after identification of consequnces of isoform switches (via \code{analyzeSwitchConsequences}) have been done as then these data are not longer needed. It is highly recomended to save a version of the switchAnalyzeRlist with all the annotation before removeing it - typically done with \code{base::saveRDS} or \code{base::save}.
}
\usage{
removeAnnoationData(
    switchAnalyzeRlist,
    removeBioSeq = TRUE,
    removeQuantData = TRUE
)
}

\arguments{
\item{switchAnalyzeRlist}{
A \code{switchAnalyzeRlist} object where \code{analyzeSwitchConsequences()} have been run to identify consequences of isoform switches
}

\item{removeBioSeq}{
A logical indicating whether to remove biological sequences (nucloetide and amino acid sequences) from the switchAnalyzeRlist. Default is TRUE.
}

\item{removeQuantData}{
A logic indicating whether to remove replicate quantification data (counts, abundances and isoform fractions) from the switchAnalyzeRlist. Default is TRUE.
}

}

\value{
This function returns a data.frame where the first column is the gene id followed by the gene expression in all samples.
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
#       "/myAnnotation/isoformsQuantified.gtf" to the functions
# 2) importRdata directly supports import of a GTF file - just supply the
#       path (e.g. "/myAnnotation/isoformsQuantified.gtf") to the isoformExonAnnoation argument

### Import quantifications
salmonQuant <- importIsoformExpression(system.file("extdata/", package="IsoformSwitchAnalyzeR"))


### Make design matrix
myDesign <- data.frame(
    sampleID = colnames(salmonQuant$abundance)[-1],
    condition = gsub('_.*', '', colnames(salmonQuant$abundance)[-1])
)

### Create switchAnalyzeRlist
aSwitchList <- importRdata(
    isoformCountMatrix   = salmonQuant$counts,
    isoformRepExpression = salmonQuant$abundance,
    designMatrix         = myDesign,
    isoformExonAnnoation = system.file("extdata/example.gtf.gz", package="IsoformSwitchAnalyzeR")
)
names(aSwitchList)

### Remove annotation data
aSwitchList <- removeAnnoationData(aSwitchList)
names(aSwitchList)


}
