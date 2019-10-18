\name{analyzeSwitchConsequences}
\alias{analyzeSwitchConsequences}
\title{
Analyze Consequences of Isoform Switches
}
\description{
This function extracts all isoforms with an absolute dIF change larger than \code{dIFcutoff} from genes with a significant isoform switch (as defined by \code{alpha}). For each gene these isoforms are then analyzed for differences in the functional annotation (defined by \code{consequencesToAnalyze}) by pairwise comparing the isoforms that are used more (switching up (dIF > 0)) with the isoforms that are used less (switching down (dIF < 0)). For each comparison a small report of the analyzed features is returned.
}
\usage{
analyzeSwitchConsequences(
    switchAnalyzeRlist,
    consequencesToAnalyze=c(
        'intron_retention',
        'coding_potential',
        'ORF_seq_similarity',
        'NMD_status',
        'domains_identified',
        'IDR_identified',
        'IDR_type',
        'signal_peptide_identified'
    ),
    alpha=0.05,
    dIFcutoff=0.1,
    onlySigIsoforms=FALSE,
    ntCutoff=50,
    ntFracCutoff=NULL,
    ntJCsimCutoff=0.8,
    AaCutoff=10,
    AaFracCutoff=0.5,
    AaJCsimCutoff=0.9,
    removeNonConseqSwitches=TRUE,
    showProgress=TRUE,
    quiet=FALSE
)
}

\arguments{
    \item{switchAnalyzeRlist}{A \code{switchAnalyzeRlist} object containing the result of an isoform switch analysis (such as the one provided by \code{isoformSwitchTestDEXSeq}) as well as additional annotation data for the isoforms.}

    \item{consequencesToAnalyze}{A vector of strings indicating what type of functional consequences to analyze. Do note that there is bound to be some differences between isoforms (else they would be identical and not annotated as separate isoforms). See details for full list of usable strings and their meaning. Default is c('intron_retention','coding_potential','ORF_seq_similarity','NMD_status','domains_identified','signal_peptide_identified') (corresponding to analyze: intron retention, CPAT result, ORF AA sequence similarity, NMD status, protein domains annotated and signal peptides annotated by Pfam).}

    \item{alpha}{The cutoff which the FDR correct p-values (q-values) must be smaller than for calling significant switches. Default is 0.05.}

    \item{dIFcutoff}{The cutoff which the changes in (absolute) isoform usage must be larger than before an isoform is considered switching. This cutoff can remove cases where isoforms with (very) low dIF values are deemed significant and thereby included in the downstream analysis. This cutoff is analogous to having a cutoff on log2 fold change in a normal differential expression analysis of genes to ensure the genes have a certain effect size. Default is 0.1 (10\%).}

    \item{onlySigIsoforms}{ A logic indicating whether to only consider significant isoforms, meaning only analyzing genes where at least two isoforms which both have significant usage changes in opposite direction (quite strict) Naturally this only works if the isoform switch test used have isoform resolution (which the build in \link{isoformSwitchTestDEXSeq} has). If FALSE all isoforms with an absolute dIF value larger than \code{dIFcutoff} in a gene with significant switches (defined by \code{alpha} and \code{dIFcutoff}) are included in the pairwise comparison. Default is FALSE (non significant isoforms are also considered based on the logic that if one isoform changes it contribution - there must be an equivalent opposite change in usage in the other isoforms from that gene).}

    \item{ntCutoff}{An integer indicating the length difference (in nt) a comparison must be larger than for reporting differences when evaluating 'isoform_length', 'ORF_length', '5_utr_length', '3_utr_length', 'isoform_seq_similarity', '5_utr_seq_similarity' and '3_utr_seq_similarity'. Default is 50 (nt).}

    \item{ntFracCutoff}{An numeric indicating the cutoff in length difference, measured as a fraction of the length of the downregulated isoform, a comparison must be larger than for reporting differences when evaluating 'isoform_length', 'ORF_length', '5_utr_length', '3_utr_length'. For example does 0.05 mean the upregulated isoform must be 5\% longer/shorter before it is reported. NULL disables the filter. Default is NULL.}

    \item{ntJCsimCutoff}{An numeric (between 0 and 1) indicating the cutoff on Jaccard Similarity (JCsim) (see details) between the overlap of two nucloetide (nt) sequences. If the measured JCsim is smaller than this cutoff the sequences are considered different and reported as such. This cutoff affects the result of the 'isoform_seq_similarity', '5_utr_seq_similarity' and '3_utr_seq_similarity' analysis. Default is 0.8}

    \item{AaCutoff}{ An integer indicating the length difference (in AA) a comparison must be larger than for reporting differences when evaluating 'ORF_seq_similarity', primarily implemented to avoid differences in very short AA sequences being classified as different. Default is 10 (AA).}

    \item{AaFracCutoff}{An numeric indicating the cutoff of length difference of the protein domain or IDR. The difference is measured as a fraction of the longest region, a comparison must be larger than before reporting it. Only used when analyzing 'domain_length' or 'IDR_length'. For example setting AaFracCutoff = 0.5 means the short protein domain (or IDR) must be <50\% of the length of the long region, before it is reported. NULL disables the filter. Default is 0.5.}

    \item{AaJCsimCutoff}{An numeric (between 0 and 1) indicating the cutoff on Jaccard Similarity (JCsim) (see details) between the overlap of two amino acid (AA) sequences. If the measured JCsim is smaller than this cutoff the sequences are considered different and reported as such. This cutoff affect the result of the 'ORF_seq_similarity' analysis. Default is 0.9}

    \item{removeNonConseqSwitches}{A logic indicating whether to remove the comparison of isoforms where no consequences were found (if TRUE) or to keep then (if FALSE). Defaults is TRUE.}

    \item{showProgress}{A logic indicating whether to make a progress bar (if TRUE) or not (if FALSE). Default is TRUE.}

    \item{quiet}{ A logic indicating whether to avoid printing progress messages (incl. progress bar). Default is FALSE}
}

\details{
Changes in isoform usage are measure as the difference in isoform fraction (dIF) values, where isoform fraction (IF) values are calculated as <isoform_exp> / <gene_exp>.\cr

The idea is that once we know there is (at least) one isoform with a significant change in how much it is used (as defined by \code{alpha} and \code{dIFcutoff}) in a gene we take that/those isoform(s) and compare the functional annotation of this isoform to the isoform(s) with the compensatory change(s) in isoform usage (since if one isoform is use more another/others have to be used less). Here we only require that one of the isoforms in the comparison of annotation is significant (unless \code{onlySigIsoforms=TRUE}, then both must be), but all isoforms considered must have a change in isoform usage larger than \code{dIFcutoff}.

Note that sometimes we find complex switches meaning that many isoforms passes all the filters. In these cases we compare all pairwise combinations of the isoform(s) used more (positive dIF) vs the isoform(s) used less (negative dIF).\cr

For sequences similarity analysis the two compared sequences are (globally) aligned to one another and the Jaccard similarity (JCsim) is calculated. Here JCsim is defined as the length of the aligned regions (omitting gaps) divided by the total combined unique sequence length: JCsim = (length of aligned region w.o gaps) / ( (length of sequence a) + (length of sequence b) - (length of aligned region w.o gaps) ). The pairwise alignment is done with \code{pairwiseAlignment{Biostrings}} as a Needleman-Wunsch global alignment which is guaranteed to find the optimal global alignment. The pairwise alignment is done with end gap penalties for the full sequences alignments ('isoform_seq_similarity' and 'ORF_seq_similarity') and without gap penalties for the alignment of sub-sequence ('5_utr_seq_similarity' and '3_utr_seq_similarity') by specifying type='global' and type='overlap' respectively.

If AA sequences were trimmed in the process of exporting the fasta files when using \code{extractSequence} the regions not analyzed in both isoforms will be ignored.

The arguments passed to \code{consequencesToAnalyze} must be a combination of:
\itemize{
    \item \code{all} : Test transcripts for any of the differences described below. Please note that jointly the analysis below covers all transcript feature meaning that they should be different. Furthermore note that 'class_code' will only be included if the switchAnalyzeRlist was made from Cufflinks/Cuffdiff output.

    \item \code{tss} : Test transcripts for whether they use different Transcription Start Site (TSS) (more than \code{ntCutoff} from each other).

    \item \code{tts} : Test transcripts for whether they use different Transcription Termination Site (TTS) (more than \code{ntCutoff} from each other).

    \item \code{last_exon} : Test whether transcripts utilizes different last exons (defined as the last exon of each transcript is non-overlapping).

    \item \code{isoform_seq_similarity} : Test whether the isoform nucleotide sequences are different (as described above). Reported as different if the measured JCsim is smaller than \code{ntJCsimCutoff} and the length difference of the aligned and combined region is larger than \code{ntCutoff}.

    \item \code{isoform_length} : Test transcripts for differences in isoform length. Only reported if the difference is larger than indicated by the \code{ntCutoff} and \code{ntFracCutoff}. Please note that this is a less powerful analysis than implemented in 'isoform_seq_similarity' as two equally long sequences might be very different.

    \item \code{exon_number} : Test transcripts for differences in exon number.

    \item \code{intron_structure} : Test transcripts for differences in intron structure, e.g. usage of exon-exon junctions. This analysis corresponds to analyzing whether all introns in one isoform is also found in the other isoforms.

    \item \code{intron_retention} : Test for differences in intron retentions (and their genomic positions). Require that \code{analyzeIntronRetention} have been run.

    \item \code{isoform_class_code} : Test transcripts for differences in the transcript classification provide by cufflinks. For a updated list of class codes see \url{http://cole-trapnell-lab.github.io/cufflinks/cuffcompare/#transfrag-class-codes}.

    \item \code{coding_potential} : Test transcripts for differences in coding potential, as indicated by the CPAT or CPC2 analysis. Requires that \code{analyzeCPAT} or \code{analyzeCPC2} have been used to add external conding potential analysis to the \code{switchAnalyzeRlist}.

    \item \code{ORF_seq_similarity} : Test whether the amino acid sequences of the ORFs are different (as described above). Reported as different if the measured JCsim is smaller than \code{AaJCsimCutoff} and the length difference of the aligned and combined region is larger than \code{AaCutoff}. Requires that least one of the isoforms are annotated with a ORF either via \code{identifyORF} or by supplying a GTF file and setting \code{addAnnotatedORFs=TRUE} when creating the switchAnalyzeRlist.

    \item \code{ORF_genomic} : Test transcripts for differences in genomic position of the Open Reading Frames (ORF). Requires that least one of the isoforms are annotated with an ORF either via \code{identifyORF} or by supplying a GTF file and setting\code{addAnnotatedORFs=TRUE} when creating the switchAnalyzeRlist.

    \item \code{ORF_length} : Test transcripts for differences in length of Open Reading Frames (ORF). Note that this is a less powerful analysis than implemented in \code{ORF_seq_similarity} as two equally long sequences might be very different. Only reported if the difference is larger than indicated by the \code{ntCutoff} and \code{ntFracCutoff}. Requires that least one of the isoforms are annotated with a ORF either via \code{identifyORF} or by supplying a GTF file and setting \code{addAnnotatedORFs=TRUE} when creating the switchAnalyzeRlist.

    \item \code{5_utr_seq_similarity} : Test whether the isoform nucleotide sequences of the 5' UnTranslated Region (UTR) are different (as described above). The 5'UTR is defined as the region from the transcript start to the ORF start. Reported as different if the measured JCsim is smaller than \code{ntJCsimCutoff} and the length difference of the aligned and combined region is larger than \code{ntCutoff}. Requires that both the isoforms are annotated with an ORF either via \code{identifyORF} or by supplying a GTF file and setting \code{addAnnotatedORFs=TRUE} when creating the switchAnalyzeRlist.

    \item \code{5_utr_length} : Test transcripts for differences in the length of the 5' UnTranslated Region (UTR). The 5'UTR is defined as the region from the transcript start to the ORF start. Note that this is a less powerful analysis than implemented in '5_utr_seq_similarity' as two equally long sequences might be very different. Only reported if the difference is larger than indicated by the \code{ntCutoff} and \code{ntFracCutoff}. Requires that both the isoforms are annotated with a ORF either via \code{identifyORF} or by supplying a GTF file and setting \code{addAnnotatedORFs=TRUE} when creating the switchAnalyzeRlist.

    \item \code{3_utr_seq_similarity} : Test whether the isoform nucleotide sequences of the 3' UnTranslated Region (UTR) are different (as described above). The 3'UTR is defined as the region from the end of the ORF to the transcript end. Reported as different if the measured JCsim is smaller than \code{ntJCsimCutoff} and the length difference of the aligned and combined region is larger than \code{ntCutoff}. Requires that both the isoforms are annotated with a ORF either via \code{identifyORF} or by supplying a GTF file and setting \code{addAnnotatedORFs=TRUE} when creating the switchAnalyzeRlist.

    \item \code{3_utr_length} : Test transcripts for differences in the length of the 3' UnTranslated Regions (UTR). The 3'UTR is defined as the region from the end of the ORF to the transcript end. Note that this is a less powerful analysis than implemented in \code{3_utr_seq_similarity} as two equally long sequences might be very different. Requires that \code{identifyORF} have been used to predict NMD sensitivity or that the ORF was imported though one of the dedicated import functions implemented in isoformSwitchAnalyzeR. Only reported if the difference is larger than indicated by the \code{ntCutoff} and \code{ntFracCutoff}. Requires that both the isoforms are annotated with a ORF either via \code{identifyORF} or by supplying a GTF file and setting \code{addAnnotatedORFs=TRUE} when creating the switchAnalyzeRlist.

    \item \code{NMD_status} : Test transcripts for differences in sensitivity to Nonsense Mediated Decay (NMD). Requires that both the isoforms have been annotated with PTC either via \code{identifyORF} or by supplying a GTF file and setting \code{addAnnotatedORFs=TRUE} when creating the switchAnalyzeRlist.

    \item \code{domains_identified} : Test transcripts for differences in the name and order of which domains are identified by the Pfam in the transcripts. Requires that \code{analyzePFAM} have been used to add external Pfam analysis to the \code{switchAnalyzeRlist}. Requires that both the isoforms are annotated with a ORF either via \code{identifyORF} or by supplying a GTF file and setting \code{addAnnotatedORFs=TRUE} when creating the switchAnalyzeRlist.

    \item \code{domain_length} : Test transcripts for differences in the length of overlapping domains of the same type (same hmm_name) thereby enabling analysis of protein domain truncation. Do however note that a small difference in length is will likely not truncate the protein domain. The length difference, measured in AA, must be larger than \code{AaCutoff} and \code{AaFracCutoff} .Requires that \code{analyzePFAM} have been used to add external Pfam analysis to the \code{switchAnalyzeRlist}. Requires that both the isoforms are annotated with a ORF either via \code{identifyORF} or by supplying a GTF file and setting \code{addAnnotatedORFs=TRUE} when creating the switchAnalyzeRlist.

    \item \code{genomic_domain_position} : Test transcripts for differences in the genomic position of the domains identified by the Pfam analysis (Will be different unless the two isoforms have the same domains at the same genomic location). Requires that \code{analyzePFAM} have been used to add external Pfam analysis to the \code{switchAnalyzeRlist}. Requires that both the isoforms are annotated with a ORF either via \code{identifyORF} or by supplying a GTF file and setting \code{addAnnotatedORFs=TRUE} when creating the switchAnalyzeRlist (and are thereby also affected by \code{removeNoncodinORFs=TRUE} in \code{analyzeCPAT}).

    \item \code{IDR_identified} : Test for differences in isoform IDRs. Specifically the two isoforms are analyzed for whether they contain IDRs which do not overlap in genomic coordinates. Requires that \code{analyzeNetSurfP2} or \code{analyzeIUPred2A} have been used to add external IDR analysis to the \code{switchAnalyzeRlist}.

    \item \code{IDR_length} : Test for differences in the length of overlapping (in genomic coordinates) IDRs. The length difference, measured in AA, must be larger than \code{AaCutoff} and \code{AaFracCutoff}. Requires that \code{analyzeNetSurfP2} or \code{analyzeIUPred2A} have been used to add external IDR analysis to the \code{switchAnalyzeRlist}.

    \item \code{IDR_type} : Test for differences in IDR type. Specifically the two isoforms are tested for overlapping IDRs (genomic coordinates) and overlapping IDRs are compared with regards to their IDR type (IDR vs IDR w binding site). Only available if \code{analyzeIUPred2A} was used to add external IDR analysis to the \code{switchAnalyzeRlist}.

    \item \code{signal_peptide_identified} : Test transcripts for differences in whether a signal peptide was identified or not by the SignalP analysis. Requires that \code{analyzeSignalP} have been used to add external SignalP analysis to the \code{switchAnalyzeRlist}. Requires that both the isoforms are annotated with a ORF either via \code{analyzeORF} or by supplying a GTF file and setting \code{addAnnotatedORFs=TRUE} when creating the switchAnalyzeRlist (and are thereby also affected by \code{removeNoncodinORFs=TRUE} in \code{analyzeCPAT}).
}
}

\value{
The supplied \code{switchAnalyzeRlist} is returned, but now annotated with the predicted functional consequences as follows. First a column called 'switchConsequencesGene' is added to \code{isoformFeatures} entry of \code{switchAnalyzeRlist}. This column containing a binary indication (TRUE/FALSE (and NA)) of whether the switching gene have predicted functional consequences or not.

Secondly the data.frame 'switchConsequence' is added to the \code{switchAnalyzeRlist} containing one row feature analyzed per comparison of isoforms pr comparison of condition. It contains 8 columns:
\itemize{
    \item{\code{gene_ref} : A unique reference to a specific gene in a specific comparison of conditions. Enables easy handles to integrate data from all the parts of a \code{switchAnalyzeRlist}.}
    \item{\code{gene_id}: The id of the gene which the isoforms compared belongs to. Matches the 'gene_id' entry in the 'isoformFeatures' entry of the switchAnalyzeRlist}
    \item{\code{gene_name} : The gene name associated with the <gene_id>, typically a more readable one (for example p53 or BRCA1) }
    \item{\code{condition_1}: The first condition of the comparison. Should be though of as the ground state - meaning the changes occure from condition_1 to condition_2. Matches the 'condition_1' entry in the 'isoformFeatures' entry of the switchAnalyzeRlist}
    \item{\code{condition_2}: The second condition of the comparison. Should be though of as the changed state - meaning the changes occure from condition_1 to condition_2. Matches the 'condition_2' entry in the 'isoformFeatures' entry of the switchAnalyzeRlist}
    \item{\code{isoformUpregulated}: The name of the isoform which is used more in condition_2 (when compared to condition_1, positive dIF values). Matches the 'isoform_id' entry in the 'isoformFeatures' entry of the switchAnalyzeRlist}
    \item{\code{isoformDownregulated}: The name of the isoform which is used less in condition_2 (when compared to condition_1, negative dIF values). Matches the 'isoform_id' entry in the 'isoformFeatures' entry of the switchAnalyzeRlist}
    \item{\code{iso_ref_up} : A unique reference to a specific isoform in a specific comparison of conditions for the isoform switching up. Enables easy handles to integrate data from all the parts of a \code{switchAnalyzeRlist}.}
    \item{\code{iso_ref_down} : A unique reference to a specific isoform in a specific comparison of conditions for the isoform switching down. Enables easy handles to integrate data from all the parts of a \code{switchAnalyzeRlist}.}
    \item{\code{featureCompared}: The category of the isoform features/annotation compared in this row (see details above)}
    \item{\code{isoformsDifferent}: A logic (TRUE/FALSE) indicating whether the two isoforms are different with respect to the featureCompared (see details above)}
    \item{\code{switchConsequence}: If the isoforms compared are different this columns contains a short description of the features of the upregulated isoform. E.g. domain loss means that the upregulated isoforms (isoformUpregulated) have lost domains compare to the downregulated isoform (isoformDownregulated).}
}
}

\references{
Vitting-Seerup et al. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. (2017).
}
\author{
Kristoffer Vitting-Seerup
}
\seealso{
\code{\link{analyzeORF}}\cr
\code{\link{analyzeCPAT}}\cr
\code{\link{analyzePFAM}}\cr
\code{\link{analyzeSignalP}}\cr
\code{\link{extractConsequenceSummary}}\cr
\code{\link{extractConsequenceEnrichment}}\cr
\code{\link{extractConsequenceEnrichmentComparison}} \cr
\code{\link{extractConsequenceGenomeWide}}
}
\examples{
### Prepare example data
data("exampleSwitchListAnalyzed")

# subset for fast runtime
exampleSwitchListAnalyzed <- subsetSwitchAnalyzeRlist(
    exampleSwitchListAnalyzed,
    exampleSwitchListAnalyzed$isoformFeatures$gene_id \%in\% sample(exampleSwitchListAnalyzed$isoformFeatures$gene_id, 10)
)

### Analyze consequences
consequencesOfInterest <- c(
    'intron_retention',
    'coding_potential',
    'NMD_status',
    'domains_identified'
)

exampleSwitchListAnalyzed <- analyzeSwitchConsequences(
    exampleSwitchListAnalyzed,
    consequencesToAnalyze = consequencesOfInterest,

)

### simple overview
extractSwitchSummary(exampleSwitchListAnalyzed, filterForConsequences = FALSE)
extractSwitchSummary(exampleSwitchListAnalyzed, filterForConsequences = TRUE)


### Detailed switch overview
consequenceSummary <- extractConsequenceSummary(
    exampleSwitchListAnalyzed,
    includeCombined = TRUE,
    returnResult = TRUE,        # return data.frame with summary
    plotGenes = TRUE            # plot summary
)



### Now switches are analyzed we can also extract the the largest/most significant switches with the extractTopSwitches() function
# Extract top 2 switching genes (by q-value)
extractTopSwitches(
    exampleSwitchListAnalyzed,
    filterForConsequences = TRUE,
    n = 2,
    extractGenes = TRUE,
    sortByQvals = TRUE
)

# Extract top 2 switching isoforms (by q-value)
extractTopSwitches(
    exampleSwitchListAnalyzed,
    filterForConsequences = TRUE,
    n = 2,
    extractGenes = FALSE,
    sortByQvals = TRUE
)

# Extract top 2 switching isoforms (by dIF)
extractTopSwitches(
    exampleSwitchListAnalyzed,
    filterForConsequences = TRUE,
    n = 2,
    extractGenes = FALSE,
    sortByQvals = FALSE
)

### Note the function ?extractConsequenceSummary is specific made for the post analysis of switching consequences
}
