\name{isoformToGeneExp}
\alias{isoformToGeneExp}
\title{
Summarize isoform expression to gene level expression.
}
\description{
Extract a data.frame with (mean) gene expression, isoform expression or Isoform Fraction values for all conditions (columns) from a switchAnalyzeRlist.
}
\usage{
isoformToGeneExp(
    isoRepExpWithGeneId,
    showProgress = TRUE,
    quiet = FALSE
)
}

\arguments{
\item{isoRepExpWithGeneId}{
A isoform expression matrix which in addition to expression columns contain two colmns 'isoform_id' and 'gene_id' indicating which isoforms are a part of which gene.
}

\item{showProgress}{
A logic indicating whether to make a progress bar (if TRUE) or not (if FALSE). Default is TRUE.
}

\item{quiet}{ A logic indicating whether to avoid printing progress messages. Default is FALSE}
}

\details{
This is just specialized function that is a lot faster than general pourpose functions such as ddplyr.
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
### Construct data needed from example data in cummeRbund package
### (The recomended way of analyzing Cufflinks/Cuffdiff datat is via importCufflinksCummeRbund
### This is jus an easy way to get some example data ).
cuffDB <- prepareCuffExample()

isoRepCount <- repCountMatrix(isoforms(cuffDB))
isoRepCount$isoform_id <- rownames(isoRepCount)

### add gene info
localAnnotaion <- import(system.file("extdata/chr1_snippet.gtf", package="cummeRbund"))[,c('transcript_id','gene_id')]
colnames(localAnnotaion@elementMetadata)[1] <- 'isoform_id'

isoRepCount$gene_id <- localAnnotaion$gene_id[match(
    isoRepCount$isoform_id , localAnnotaion$isoform_id
)]

### Summarize to gene level
geneRepCount <- isoformToGeneExp(isoRepCount)


}
