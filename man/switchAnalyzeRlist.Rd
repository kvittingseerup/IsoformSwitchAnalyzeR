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
    conditions,
    sourceId
)
}

\arguments{
	\item{isoformFeatures}{A \code{data.frame} where each row corresponds to a isoform in a specific comparison and contains all the annotation for this isoform. See \code{details} below for details.}

	\item{exons}{A \code{GRanges} object containing isoform exon structure. See \code{details} below for details.}

	\item{conditions}{A data.frame giving the number of biological replicates in each condition. See \code{details} below for details.}

	\item{sourceId}{A character stating the origin of the data used to create the switchAnalyzeRlist.}
}


\details{
For cufflinks data, use \link{importCufflinksCummeRbund} or \link{importCufflinksFiles} to prepare the switchAnalyzeRlist. For data from StringTie or RSEM use \link{importBallgownData}.
For other RNA-seq assemblies, either uses this constructor or the general-purpose \link{importRdata} to create the switchAnalyzeRlist.

\cr
The \code{isoformFeatures} should be a data.frame where each row corresponds to a isoform in a specific comparison and contains all the annoation for this isoform. The data.frame can contain any colums supplied (enabling addition of user specified columns) but the following columns are nessesary and must be provided:
    \itemize{\code{isoform_id} : A unique isoform id}
    \itemize{\code{gene_id} : A unique gene id referring to a gene at a specific genomic loci (not the same as gene_name since gene_names can refer to multiple genomic loci)}
    \itemize{\code{condition_1} : Name of the first condition  in the comparison}
    \itemize{\code{condition_2} : Name of the second condition  in the comparison}
    \itemize{\code{gene_name} : The gene name associated with the <gene_id>, typically a more readable one (for example p53 or BRCA1) }

    \itemize{\code{gene_value_1} : Expression of <gene_id> in condition_1 (must be normalized) }
    \itemize{\code{gene_value_2} : Expression of <gene_id> in condition_2 (must be normalized) }
    \itemize{\code{gene_stderr_1} : Standard error (of mean) of <gene_id> expression in condition_1 }
    \itemize{\code{gene_stderr_2} : Standard error (of mean) of <gene_id> expression in condition_2 }
    \itemize{\code{gene_log2_fold_change} : log2 fold change of <gene_id> expression between condition_1 and condition_2 }
    \itemize{\code{gene_q_value} : The FDR corrected (for multiple testing) p-value of the differential expression test of <gene_id> }

    \itemize{\code{iso_value_1} : Expression of <isoform_id> in condition_1 (must be normalized)}
    \itemize{\code{iso_value_2} : Expression of <isoform_id> in condition_2 (must be normalized)}
    \itemize{\code{iso_stderr_1} : Standard error (of mean) of <isoform_id> expression in condition_1}
    \itemize{\code{iso_stderr_2} : Standard error (of mean) of <isoform_id> expression in condition_2}
    \itemize{\code{iso_log2_fold_change} : log2 fold change of <isoform_id> expression between condition_1 and condition_2}
    \itemize{\code{iso_q_value} : The FDR corrected (for multiple testing) p-value of the differential expression test of <isoform_id>}

    \itemize{\code{IF1} : The <isoform_id> usage in condtion 1 (given as Isoform Fraction (IF) value)}
    \itemize{\code{IF2} : The <isoform_id> usage in condtion 2 (given as Isoform Fraction (IF) value)}
    \itemize{\code{dIF} : The change in isoform usage from condtion_1 to condition_2 (difference in IF values (dIF))}

    \itemize{\code{isoform_switch_q_value} : The q-value of the test of differential isoform usage in <isoform_id> between condtion 1 and condition 2. Use NA if not performed. Will be overwritten by the result of \code{testIsoformSwitches}. If only performed at gene level use same values on isoform level.}
    \itemize{\code{gene_switch_q_value} : The q-value of the test of differential isoform usage in <gene_id> between condtion 1 and condition 2. Use NA if not performed. Will be overwritten by the result of \code{testIsoformSwitches}.}

\cr
The \code{exons} argument must be supplied with a \code{GenomicRange} object containing one entry pr exon in each isoform. Furthermore it must also have two meta collumns called \code{isoform_id} and \code{gene_id} which links it to the information in the \code{isoformFeatures} entry.

\cr
The \code{conditions} should be a data.frame with two collumns: \code{condition} and \code{nrReplicates} giving the number of biological (not technical) replicates each condition analyzed. The strings used to conditions the condtions must match the strings used in \code{condition_1} and \code{condition_2} columns of the \code{isoformFeatures} entry.
}

\value{
A list-type object \code{switchAnalyzeRlist} object containing all the information needed to do the full analysis with IsoformSwitchAnalyzeR. Note that \code{switchAnalyzeRlist} appears as a normal list and all the information (incl that added by all the analyze* functions) can be obtained using both the named entries (f.x. myIsoSwitchList$isoformFeatures ) or indexes (f.x myIsoSwitchList[[1]] ).

When fully analyzed the \code{isoformFeatures} entry of the  will furthermore contain the following collumns:
\itemize{\code{id}: During the creation of \code{switchAnalyzeRlist} a unique id is constructed for each row - meaning for each isoform in each comparison. The id is constructed as 'soComp' an acronym for 'isoform comparison', followed by XXXXXXXX indicating the row number}
\itemize{\code{PTC}: A logic indicating whether the <isoform_id> is classified as having a Premature Termination Codon. This is defined as having a stopcodon more than \code{PTCDistance}(default is 50) nt upstream of the last exon exon.}
\itemize{\code{codingPotentialValue}: containing the coding potential value predicted by CPAT.}
\itemize{\code{codingPotential}: A logic (TRUE/FALSE) indicating whether the isoform is coding or not (based on the \code{codingCutoff} supplied)}
\itemize{\code{signal_peptide_identified}: A string ('yes'/'no') indicating whether the <isoform_id> have a signal peptide, as predicted by SignalP.}
\itemize{\code{domain_identified}: A string ('yes'/'no') indicating whether the <isoform_id> contain (at least one) protein domain, as predicted by pfam.}
\itemize{\code{switchConsequencesGene}: A logic (TRUE/FALSE) indicating whether the <gene_id> contain an isoform switch with functional consequences, as predicted by \code{analyzeSwitchConsequences}.}

}

\seealso{
\code{\link{importRdata}}
\code{\link{readCufflinks}}
\code{\link{importCufflinksCummeRbund}}
\code{\link{importCufflinksFiles}}
\code{\link{importGTF}}
\code{\link{importBallgownData}}
}
\author{Kristoffer Vitting-Seerup}
\references{ Manuscript in preperation}

\examples{
### A switchAnalyzeRlist needs 3 annotation types - let's create them one at the time.

# 1) Make "isoformFeatures" entry
isoAnnotation <- data.frame(
    isoform_id = paste('iso',1:3,sep='_'),
    gene_id='gene_1',
    condition_1="groundState",
    condition_2="modifiedState",
    gene_name='favoriteGene',
    gene_value_1=10,
    gene_value_2=11,
    gene_stderr_1=2,
    gene_stderr_2=2,
    gene_log2_fold_change=log2(11/10),
    gene_q_value=1,
    iso_value_1=c(1,8,1),
    iso_value_2=c(1,2,7),
    iso_stderr_1=rep(0.1,3),
    iso_stderr_2=rep(0.2,3),
    iso_log2_fold_change=log2( c(1,2,7) / c(1,8,1) ),
    iso_p_value=1,
    iso_q_value=1,
    IF1=c(1,8,1) / 10,
    IF2=c(1,2,7) / 11,
    dIF=(c(1,2,7) / 11) - (c(1,8,1) / 10),
    isoform_switch_q_value=NA,
    gene_switch_q_value=NA,
    stringsAsFactors = FALSE
)

# 2) make "conditions" entry
repOverview <- data.frame(
    condition=c('groundState','modifiedState'),
    nrReplicates=c(3,3),
    row.names = NULL, stringsAsFactors = FALSE
)

# 3) make "exons" entry
myExons <- GRanges(
    seqnames = 'chr1',
    ranges = IRanges(start = c(1,1,10), end = c(20,20,20)),
    strand = '+',
    isoform_id=paste('iso',1:3,sep='_'),
    gene_id='gene_1'
)

### Combine it all to make a switchAnalyzeRlist

aSwitchList <- createSwitchAnalyzeRlist(
    isoformFeatures=isoAnnotation,
    exons=myExons,
    conditions=repOverview,
    sourceId='homeMade'
)
aSwitchList
}
