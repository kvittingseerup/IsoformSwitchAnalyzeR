\name{SpliceRList-class}
\docType{class}
\alias{SpliceRList-class}

\title{Class \code{"SpliceRList"}}
\description{
Creats a SpliceRList - the backbone object of the R package spliceR - containing all the information nessesary for analysis performed with spliceR.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SpliceRList", <data_list>)} where <data_list> must be a list with the following named entries:

\itemize{
	\item{transcript_features}{\code{GRanges} object containing transcript features including the meta columns \code{gene_id} for gene id, \code{isoform_id} for transcript id, \code{sample_1} and \code{sample_2} for sample identifiers, \code{expression_1} and \code{expression_2} for expression values for sample 1 and sample 2, respectively (typically FPKM values or some other normalized count values)}
	\item{exon_features}{\code{GRanges} object containing transcript exons including the meta columns \code{gene_id} for gene id and \code{isoform_id} for transcript id.}
	\item{assembly_id}{character, giving genome assemlby ('hg19', 'hg18', 'mm9', etc.).}
	\item{source_id}{A character, either "cufflinks" or "granges", stating source of transcript assembly.}
	\item{conditions}{A character vector, giving the samples or conditions for the RNA-seq experiment.}
	\item{transcripts_plot}{A dataframe, reserved for plotting functions. Initialize as NULL}
	\item{filter_params}{A character vector, reserved for plotting functions. Initialize as NULL}
}
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
}
\section{Methods}{
No methods defined with class "SpliceRList" in the signature.
}
\references{
Vitting-Seerup K , Porse BT, Sandelin A, Waage J. (2014) spliceR: an R package for classification of alternative splicing and prediction of coding potential from RNA-seq data. BMC Bioinformatics 15:81.
}
\author{
Kristoffer Vitting-Seerup
}
\examples{
showClass("SpliceRList")
}
\keyword{classes}
