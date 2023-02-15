\name{createSwitchAnalyzeRlist}
\alias{createSwitchAnalyzeRlist}
\alias{switchAnalyzeRlist-class}
\alias{switchAnalyzeRlist}
\docType{class}
\title{Create a switchAnalyzeRlist Object}

\description{
Create a \code{switchAnalyzeRlist} containing all the information needed to do the full analysis with IsoformSwitchAnalyzeR.
}

\usage{
createSwitchAnalyzeRlist(
    isoformFeatures,
    exons,
    designMatrix,
    isoformCountMatrix=NULL,
    isoformRepExpression=NULL,
    sourceId,
    removeFusionTranscripts = TRUE
)
}

\arguments{
	\item{isoformFeatures}{A \code{data.frame} where each row corresponds to a isoform in a specific comparison and contains all the annotation for this isoform. See \code{details} below for details.}
	\item{exons}{A \code{GRanges} object containing isoform exon structure. See \code{details} below for details.}
	\item{designMatrix}{A data.frame with the information of which samples originate from which conditions. A data.frame with two columns: \code{sampleID 1} contains the sample names which matches the column names used in \code{isoformCountMatrix}. \code{condition}: which indicates which conditions the sample originate from. If sample 1-3 originate form the same condition they should all have the same string (for example 'ctrl', in this column). By adding additional columns to this designMatrix batch effects can be taking into account with the downstream isoform switch test.}

	\item{isoformCountMatrix}{ A data.frame with unfiltered biological (not technical) replicate isoform (estimated) counts. Must have a column called 'isoform_id' with the isoform_id that matches \code{isoformFeatures}. The name of the columns must match the sample names in the \code{designMatrix} argument and contain the estimated counts.}

	\item{isoformRepExpression}{ A data.frame with unfiltered biological (not technical) replicate isoform abundances. Must have a column called 'isoform_id' with the isoform_id that matches \code{isoformFeatures}. The name of the columns must match the sample names in the \code{designMatrix} argument and contain the estimated abundances.}

	\item{sourceId}{A character stating the origin of the data used to create the switchAnalyzeRlist.}

	\item{removeFusionTranscripts}{A logic indicating whether to remove genes with cross-chromosome fusion transcripts as IsoformSwitchAnalyzeR cannot handle them.}
}


\details{
For cufflinks data, use \link{importCufflinksFiles} to prepare the switchAnalyzeRlist.
For other RNA-seq assemblies, either uses this constructor or the general-purpose \link{importRdata} to create the switchAnalyzeRlist - See vignette for details.


The \code{isoformFeatures} should be a data.frame where each row corresponds to a isoform in a specific comparison and contains all the annotation for this isoform. The data.frame can contain any columns supplied (enabling addition of user specified columns) but the following columns are necessary and must be provided:
\itemize{
    \item{\code{iso_ref} : A unique reference to a specific isoform in a specific comparison of conditions. Mainly created to have an easy handle to integrate data from all the parts of a \code{switchAnalyzeRlist}.}
    \item{\code{gene_ref} : A unique reference to a specific gene in a specific comparison of conditions. Mainly created to have an easy handle to integrate data from all the parts of a \code{switchAnalyzeRlist}.}
    \item{\code{isoform_id} : A unique isoform id}
    \item{\code{gene_id} : A unique gene id referring to a gene at a specific genomic loci (not the same as gene_name since gene_names can refer to multiple genomic loci)}
    \item{\code{condition_1} : Name of the first condition  in the comparison}
    \item{\code{condition_2} : Name of the second condition  in the comparison}
    \item{\code{gene_name} : The gene name associated with the <gene_id>, typically a more readable one (for example p53 or BRCA1) }

    \item{\code{gene_overall_mean} : Mean expression of <gene_id> across all samples (if you create it yourself consider inter-library normalization) }
    \item{\code{gene_value_1} : Expression of <gene_id> in condition_1 (if you create it yourself consider inter-library normalization) }
    \item{\code{gene_value_2} : Expression of <gene_id> in condition_2 (if you create it yourself consider inter-library normalization) }
    \item{\code{gene_stderr_1} : Standard error (of mean) of <gene_id> expression in condition_1 }
    \item{\code{gene_stderr_2} : Standard error (of mean) of <gene_id> expression in condition_2 }
    \item{\code{gene_log2_fold_change} : log2 fold change of <gene_id> expression between condition_1 and condition_2 }
    \item{\code{gene_q_value} : The FDR corrected (for multiple testing) p-value of the differential expression test of <gene_id> }

    \item{\code{iso_overall_mean} : Mean expression of <isoform_id> across all samples (if you create it yourself consider inter-library normalization) }
    \item{\code{iso_value_1} : Expression of <isoform_id> in condition_1 (if you create it yourself consider inter-library normalization)}
    \item{\code{iso_value_2} : Expression of <isoform_id> in condition_2 (if you create it yourself consider inter-library normalization)}
    \item{\code{iso_stderr_1} : Standard error (of mean) of <isoform_id> expression in condition_1}
    \item{\code{iso_stderr_2} : Standard error (of mean) of <isoform_id> expression in condition_2}
    \item{\code{iso_log2_fold_change} : log2 fold change of <isoform_id> expression between condition_1 and condition_2}
    \item{\code{iso_q_value} : The FDR corrected (for multiple testing) p-value of the differential expression test of <isoform_id>}

    \item{\code{IF_overall} : The average <isoform_id> usage across all samples (given as Isoform Fraction (IF) value)}
    \item{\code{IF1} : The <isoform_id> usage in condition 1 (given as Isoform Fraction (IF) value)}
    \item{\code{IF2} : The <isoform_id> usage in condition 2 (given as Isoform Fraction (IF) value)}
    \item{\code{dIF} : The change in isoform usage from condition_1 to condition_2 (difference in IF values (dIF))}

    \item{\code{isoform_switch_q_value} : The q-value of the test of differential isoform usage in <isoform_id> between condition 1 and condition 2. Use NA if not performed. Will be overwritten by the result of \code{testIsoformSwitches}. If only performed at gene level use same values on isoform level.}
    \item{\code{gene_switch_q_value} : The q-value of the test of differential isoform usage in <gene_id> between condition 1 and condition 2. Use NA if not performed. Will be overwritten by the result of \code{testIsoformSwitches}.}
}


The \code{exons} argument must be supplied with a \code{GenomicRange} object containing one entry pr exon in each isoform. Furthermore it must also have two meta columns called \code{isoform_id} and \code{gene_id} which links it to the information in the \code{isoformFeatures} entry.


The \code{conditions} should be a data.frame with two columns: \code{condition} and \code{nrReplicates} giving the number of biological (not technical) replicates each condition analyzed. The strings used to conditions the conditions must match the strings used in \code{condition_1} and \code{condition_2} columns of the \code{isoformFeatures} entry.
}

\value{
A list-type object \code{switchAnalyzeRlist} object containing all the information needed to do the full analysis with IsoformSwitchAnalyzeR. Note that \code{switchAnalyzeRlist} appears as a normal list and all the information (incl that added by all the analyze* functions) can be obtained using both the named entries (f.x. myIsoSwitchList$isoformFeatures ) or indexes (f.x myIsoSwitchList[[1]] ).

When fully analyzed the \code{isoformFeatures} entry of the  will furthermore contain the following columns:
\itemize{
    \item{\code{id}: During the creation of \code{switchAnalyzeRlist} a unique id is constructed for each row - meaning for each isoform in each comparison. The id is constructed as 'isoComp' an acronym for 'isoform comparison', followed by XXXXXXXX indicating the row number}
    \item{\code{PTC}: A logic indicating whether the <isoform_id> is classified as having a Premature Termination Codon. This is defined as having a stopcodon more than \code{PTCDistance}(default is 50) nt upstream of the last exon exon.}
    \item{\code{codingPotentialValue}: containing the coding potential value predicted by CPAT.}
    \item{\code{codingPotential}: A logic (TRUE/FALSE) indicating whether the isoform is coding or not (based on the \code{codingCutoff} supplied)}
    \item{\code{signal_peptide_identified}: A string ('yes'/'no') indicating whether the <isoform_id> have a signal peptide, as predicted by SignalP.}
    \item{\code{domain_identified}: A string ('yes'/'no') indicating whether the <isoform_id> contain (at least one) protein domain, as predicted by pfam.}
    \item{\code{switchConsequencesGene}: A logic (TRUE/FALSE) indicating whether the <gene_id> contain an isoform switch with functional consequences, as predicted by \code{analyzeSwitchConsequences}.}
}
}

\seealso{
\code{\link{importRdata}}\cr
\code{\link{importCufflinksFiles}}\cr
\code{\link{importGTF}}\cr
\code{\link{importIsoformExpression}}
}
\author{Kristoffer Vitting-Seerup}
\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
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
aSwitchList
}
