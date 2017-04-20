\name{SpliceRList}
\alias{SpliceRList}
\alias{SpliceRList-class}
\docType{class}
\title{
	Transcript data and annotation object for spliceR
}
\description{
Creates a \code{SpliceRList} object from two \code{GRanges} objects, an assembly id, and a source id.
The first \code{GRanges}, transcript_features, containing a list of transcripts, and including the columns \code{gene_id} for gene id, \code{tx_id} for transcript id, \code{sample_1} and \code{sample_2} for sample identifiers, \code{expression_1} and \code{expression_2} for expression values for sample 1 and sample 2, respectively (typically FPKM values or some other normalized count values), and additional optional columns (see \code{prepareCuff}).
The second, exon_features, containing a list of exons, and including the columns \code{gene_id} for gene id and \code{tx_id} for transcript id.
Assembly id, denoting genome assembly ('hg19', 'hg18', 'mm9', etc.)
Source id, denoting source of transcript assembly (currently 'cufflinks' or 'other')
Note, that the cromosome identifiers should match the assembly. For experiments 

}
\usage{
SpliceRList(transcript_features, exon_features, assembly_id, source_id, conditions, transcripts_plot=NULL,filter_params=NULL)
}
\arguments{
	\item{transcript_features}{\code{GRanges} object containing transcript features.}
	\item{exon_features}{\code{GRanges} object containing transcript features.}
	\item{assembly_id}{character, giving genome assemlby.}
	\item{source_id}{A character, either "cufflinks" or "granges", stating source of transcript assembly.}
	\item{conditions}{A character vector, giving the samples or conditions for the RNA-seq experiment.}
	\item{transcripts_plot}{A dataframe, reserved for plotting functions}
	\item{filter_params}{A character vector, reserved for plotting functions.}
}
\details{
For cufflinks data, call \code{prepareCuff} to prepare a SpliceRList.
For other RNA-seq assemblies, use this constructor to create a SpliceRList.

See the spliceR vignette for an example of creating a spliceRList from another source than Cufflinks.
}

\value{
A \code{SpliceRList} object.
}

\author{Kristoffer Knudsen, Johannes Waage}

\references{
Vitting-Seerup K , Porse BT, Sandelin A, Waage J. (2014) spliceR: an R package for classification of alternative splicing and prediction of coding potential from RNA-seq data. BMC Bioinformatics 15:81.
}


\keyword{classes}
