### For analyzing consequences
analyzeSwitchConsequences <- function(
    switchAnalyzeRlist,
    consequencesToAnalyze = c(
        'intron_retention',
        'coding_potential',
        'ORF_seq_similarity',
        'NMD_status',
        'domains_identified',
        'IDR_identified',
        'IDR_type',
        'signal_peptide_identified'
    ),
    alpha = 0.05,
    dIFcutoff = 0.1,
    onlySigIsoforms = FALSE,
    ntCutoff = 50,
    ntFracCutoff = NULL,
    ntJCsimCutoff = 0.8,
    AaCutoff = 10,
    AaFracCutoff = 0.5,
    AaJCsimCutoff = 0.9,
    removeNonConseqSwitches = TRUE,
    showProgress = TRUE,
    quiet = FALSE
) {
    ### Check input
    if (TRUE) {
        # check switchAnalyzeRlist
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(
                'The object supplied to \'switchAnalyzeRlist\' is not a \'switchAnalyzeRlist\''
            )
        }

        if (alpha < 0 |
            alpha > 1) {
            warning('The alpha parameter should usually be between 0 and 1 ([0,1]).')
        }
        if (alpha > 0.05) {
            warning(
                'Most journals and scientists consider an alpha larger than 0.05 untrustworthy. We therefore recommend using alpha values smaller than or queal to 0.05'
            )
        }

        # test wether switching have been analyzed
        if (!any(!is.na(
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value
        ))) {
            stop(
                'The analsis of isoform switching must be performed before functional consequences can be analyzed. Please run ?detectIsoformSwitching and try again.'
            )
        }

        acceptedTypes <- c(
            # Transcript
            'tss',
            'tts',
            'last_exon',
            'isoform_length',
            'exon_number',
            'intron_structure',
            'intron_retention',
            'isoform_class_code',
            # cpat
            'coding_potential',
            # ORF
            'ORF_genomic',
            'ORF_length',
            '5_utr_length',
            '3_utr_length',
            # seq similarity
            'isoform_seq_similarity',
            'ORF_seq_similarity',
            '5_utr_seq_similarity',
            '3_utr_seq_similarity',
            # ORF
            'NMD_status',
            # pfam
            'domains_identified',
            'genomic_domain_position',
            'domain_length',

            # SignalIP
            'signal_peptide_identified',

            # IDR
            'IDR_identified',
            'IDR_length',
            'IDR_type',

            # sub cell
            'sub_cell_location',
            'solubility_status'
        )

        if (!all(consequencesToAnalyze %in% c('all', acceptedTypes))) {
            stop(
                paste(
                    'The argument(s) supplied to \'typeOfconsequence\' are not accepted.',
                    'Please see ?analyzeSwitchConsequences for description of which strings are allowed.',
                    'The problem is:',
                    paste(setdiff(
                        consequencesToAnalyze , c('all', acceptedTypes)
                    ), collapse = ', '),
                    sep = ' '
                )
            )
        }

        if ('all' %in% consequencesToAnalyze) {
            consequencesToAnalyze <- acceptedTypes
        }

        ## Test whether annotation is advailable
        if ('intron_retention'  %in% consequencesToAnalyze) {
            if (is.null(switchAnalyzeRlist$intronRetentionAnalysis) & is.null( switchAnalyzeRlist$AlternativeSplicingAnalysis)) {
                stop(
                    'To test for intron retention alternative splicing must first be classified. Please run analyzeIntronRetention() and try again.'
                )
            }
        }

        if (grepl('cufflinks', switchAnalyzeRlist$sourceId)) {
            if ('isoform_class_code'  %in% consequencesToAnalyze) {
                if (!'class_code' %in%
                    colnames(switchAnalyzeRlist$isoformFeatures)
                ) {
                    stop(
                        'The switchAnalyzeRlist does not contail the calss_code information'
                    )
                }
            }
        } else {
            consequencesToAnalyze <-
                setdiff(consequencesToAnalyze, 'isoform_class_code')
        }

        if (any(
            c('ORF_genomic', 'ORF_length', 'NMD_status')  %in%
            consequencesToAnalyze
        )) {
            if (!'PTC' %in% colnames(switchAnalyzeRlist$isoformFeatures)) {
                stop(
                    'To test differences in ORFs or PCT, ORF must be annotated. Please run annotateORF() and try again'
                )
            }
        }
        if ('coding_potential'  %in% consequencesToAnalyze) {
            if (
                !'codingPotential' %in%
                colnames(switchAnalyzeRlist$isoformFeatures))
            {
                stop(
                    'To test differences in coding_potential, the result of the CPAT analysis must be advailable. Please run analyzeCPAT() or analyzeCPC2 and try again.'
                )
            }
        }
        if (any(
            c(
                'domains_identified',
                'domain_length',
                'genomic_domain_position'
            )  %in% consequencesToAnalyze
        )) {
            if (is.null(switchAnalyzeRlist$domainAnalysis)) {
                stop(
                    'To test differences in protein domains, the result of the Pfam analysis must be advailable. Please run analyzePFAM() and try again.'
                )
            }
        }
        if ('signal_peptide_identified'  %in% consequencesToAnalyze) {
            if (is.null(switchAnalyzeRlist$signalPeptideAnalysis)) {
                stop(
                    'To test differences in signal peptides, the result of the SignalP analysis must be advailable. Please run analyzeSignalP() and try again.'
                )
            }
        }
        if ( any(c('IDR_identified','IDR_type')  %in% consequencesToAnalyze)) {
            if (is.null(switchAnalyzeRlist$idrAnalysis)) {
                stop(
                    'To test differences in IDR, the result of the NetSurfP2 analysis must be advailable. Please run analyzeNetSurfP2() and try again,'
                )
            }
        }
        if( 'IDR_type' %in% consequencesToAnalyze ) {
            if( ! 'idr_type' %in% colnames(switchAnalyzeRlist$idrAnalysis) ) {
                stop('To analyse IDR_type the IDR analysis must have been done using IUPred2A and imported with the analyzeIUPred2A() function.')
            }
        }

        if (!is.numeric(ntCutoff)) {
            stop('The \'ntCutoff\' arugment must be an numeric')
        }
        if (ntCutoff <= 0) {
            stop('The \'ntCutoff\' arugment must be an numeric > 0')
        }

        if (!is.null(ntFracCutoff)) {
            if (ntFracCutoff <= 0 | ntFracCutoff > 1) {
                stop(
                    'The \'ntFracCutoff\' arugment must be a numeric in the interval (0,1]. Use NULL to disable.'
                )
            }
        }

        ### test sequence annotation
        if (any(
            consequencesToAnalyze %in% c(
                'isoform_seq_similarity',
                '5_utr_seq_similarity',
                '3_utr_seq_similarity'
            )
        )) {
            if (!any(names(switchAnalyzeRlist) == 'ntSequence')) {
                stop(
                    'The transcrip nucleotide sequences must be added to the switchAnalyzeRlist before overlap analysis can be performed. Please run \'extractSequence\' and try again.'
                )
            }
        }
        if (any(consequencesToAnalyze %in% c('ORF_seq_similarity'))) {
            if (!any(names(switchAnalyzeRlist) == 'aaSequence')) {
                stop(
                    'The transcrip ORF amino acid sequences must be added to the switchAnalyzeRlist before ORF overlap analysis can be performed. Please run \'extractSequence\' and try again.'
                )
            }
        }

        if ('sub_cell_location'  %in% consequencesToAnalyze) {
            if ( ! 'sub_cell_location' %in% colnames(switchAnalyzeRlist$isoformFeatures)) {
                stop(
                    'Cannot test for differences in sub-cellular location as such results are not annotated'
                )
            }
        }
        if ('solubility_status'  %in% consequencesToAnalyze) {
            if ( ! 'solubility_status' %in% colnames(switchAnalyzeRlist$isoformFeatures)) {
                stop(
                    'Cannot test for differences in solubility status as such results are not annotated'
                )
            }
        }

    }

    if (showProgress & !quiet) {
        progressBar <- 'text'
    } else {
        progressBar <- 'none'
    }

    ### Subset to relevant data
    if (TRUE) {
        if (!quiet) {
            message('Step 1 of 4: Extracting genes with isoform switches...')
        }

        ### Extract Iso pairs
        pairwiseIsoComparison <- extractSwitchPairs(
            switchAnalyzeRlist,
            alpha = alpha,
            dIFcutoff = dIFcutoff,
            onlySigIsoforms = onlySigIsoforms
        )

        ### Extract isoform_id pairs
        pairwiseIsoComparisonUniq <-
            unique(pairwiseIsoComparison[,c(
                'isoformUpregulated', 'isoformDownregulated'
            )])
        pairwiseIsoComparisonUniq$comparison <-
            1:nrow(pairwiseIsoComparisonUniq)

        ### Generate size reduced switchAnalyzeRList
        minimumSwitchList <- makeMinimumSwitchList(
            orgSwitchList = switchAnalyzeRlist,
            isoformsToKeep = unique(
                c(
                    pairwiseIsoComparisonUniq$isoformUpregulated,
                    pairwiseIsoComparisonUniq$isoformDownregulated
                )
            ))


    }

    ### Loop over all the the resulting genes and do a all pairwise comparison between up and down.
    if (TRUE) {
        if (!quiet) {
            message(
                paste(
                    'Step 2 of 4: Analyzing',
                    nrow(pairwiseIsoComparisonUniq),
                    'pairwise isoforms comparisons...',
                    sep = ' '
                )
            )
        }

        consequencesOfIsoformSwitching <- plyr::dlply(
            .data = pairwiseIsoComparisonUniq,
            .variables = 'comparison',
            .parallel = FALSE,
            .inform = TRUE,
            .progress = progressBar,
            .fun = function(aDF) {
                # aDF <- pairwiseIsoComparisonUniq[246,]
                compareAnnotationOfTwoIsoforms(
                    switchAnalyzeRlist    = minimumSwitchList,
                    consequencesToAnalyze = consequencesToAnalyze,
                    upIso                 = aDF$isoformUpregulated,
                    downIso               = aDF$isoformDownregulated,
                    ntCutoff              = ntCutoff,
                    ntFracCutoff          = ntFracCutoff,
                    ntJCsimCutoff         = ntJCsimCutoff,
                    AaCutoff              = AaCutoff,
                    AaFracCutoff          = AaFracCutoff,
                    AaJCsimCutoff         = AaJCsimCutoff,
                    addDescription        = TRUE,
                    testInput             = FALSE # already done by this function
                )
            }
        )

        ### Remove to those instances where there where no consequences
        if (removeNonConseqSwitches) {
            ## Remove to those instances where there where no consequences
            consequencesOfIsoformSwitching <-
                consequencesOfIsoformSwitching[which(
                    sapply(consequencesOfIsoformSwitching, function(aDF) {
                        any(aDF$isoformsDifferent)
                    })
                )]
            if (!length(consequencesOfIsoformSwitching)) {
                stop('No isoform switches with the analyzed consequences were found.')
            }
        }
    }

    ### Massage result
    if (TRUE) {
        if (!quiet) {
            message(paste(
                'Step 3 of 4: Massaging isoforms comparisons results...',
                sep = ' '
            ))
        }
        ### Convert from list to df
        consequencesOfIsoformSwitchingDf <-
            myListToDf(consequencesOfIsoformSwitching)

        ### Add the origin info
        consequencesOfIsoformSwitchingDfcomplete <- dplyr::inner_join(
            consequencesOfIsoformSwitchingDf,
            pairwiseIsoComparison,
            by = c('isoformUpregulated', 'isoformDownregulated')
        )
        #### Add additional information
        #consequencesOfIsoformSwitchingDfcomplete <- dplyr::inner_join(
        #    consequencesOfIsoformSwitchingDfcomplete,
        #    switchAnalyzeRlist$isoformFeatures[match(
        #        unique(consequencesOfIsoformSwitchingDfcomplete$gene_ref),
        #        switchAnalyzeRlist$isoformFeatures$gene_ref
        #    ),
        #    c('gene_ref',
        #      'gene_id',
        #      'gene_name',
        #      'condition_1',
        #      'condition_2')],
        #    by = 'gene_ref'
        #)

        ### reorder
        newOrder <- na.omit(match(
            c(
                'gene_ref',
                'gene_id',
                'gene_name',
                'condition_1',
                'condition_2',
                'isoformUpregulated',
                'isoformDownregulated',
                'iso_ref_up',
                'iso_ref_down',
                'featureCompared',
                'isoformsDifferent',
                'switchConsequence'
            ),
            colnames(consequencesOfIsoformSwitchingDfcomplete)
        ))
        consequencesOfIsoformSwitchingDfcomplete <-
            consequencesOfIsoformSwitchingDfcomplete[, newOrder]

        consequencesOfIsoformSwitchingDfcomplete <-
            consequencesOfIsoformSwitchingDfcomplete[order(
                consequencesOfIsoformSwitchingDfcomplete$gene_ref,
                consequencesOfIsoformSwitchingDfcomplete$isoformUpregulated,
                consequencesOfIsoformSwitchingDfcomplete$isoformDownregulated
            ), ]

    }

    ### Add result to switchAnalyzeRlist
    if (TRUE) {
        if (!quiet) {
            message('Step 4 of 4: Preparing output...')
        }
        ### Add full analysis
        switchAnalyzeRlist$switchConsequence <-
            consequencesOfIsoformSwitchingDfcomplete

        # extract indexes of those analyzed
        indexesAnalyzed <-
            which(
                switchAnalyzeRlist$isoformFeatures$gene_ref %in%
                    pairwiseIsoComparison$gene_ref
            )

        # Add indicator of switch consequence to those analyzed
        switchAnalyzeRlist$isoformFeatures$switchConsequencesGene <- NA

        switchAnalyzeRlist$isoformFeatures$switchConsequencesGene[
            indexesAnalyzed
        ] <- switchAnalyzeRlist$isoformFeatures$gene_ref[indexesAnalyzed] %in%
            consequencesOfIsoformSwitchingDfcomplete$gene_ref[which(
                consequencesOfIsoformSwitchingDfcomplete$isoformsDifferent
            )]
    }

    if (!quiet) {
        totalNrGene <-
            extractSwitchSummary(
                switchAnalyzeRlist,
                filterForConsequences = TRUE
            )
        totalNrGene <- totalNrGene$nrGenes[which(
            totalNrGene$Comparison == 'combined'
        )]
        message(
            paste(
                'Identified',
                totalNrGene,
                'genes with containing isoforms switching with functional consequences...',
                sep = ' '
            )
        )
    }
    return(switchAnalyzeRlist)
}

compareAnnotationOfTwoIsoforms <- function(
    switchAnalyzeRlist,
    consequencesToAnalyze = 'all',
    upIso,
    downIso,
    addDescription = TRUE,
    onlyRepportDifferent = FALSE,
    ntCutoff = 50,
    ntFracCutoff = NULL,
    ntJCsimCutoff = 0.8,
    AaCutoff = 10,
    AaFracCutoff = 0.5,
    AaJCsimCutoff = 0.9,
    testInput = TRUE
) {
    ### Test and prepare input
    if (testInput) {
        # check switchAnalyzeRlist
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(
                'The object supplied to \'switchAnalyzeRlist\' is not a \'switchAnalyzeRlist\''
            )
        }

        if (alpha < 0 |
            alpha > 1) {
            warning('The alpha parameter should usually be between 0 and 1 ([0,1]).')
        }
        if (alpha > 0.05) {
            warning(
                'Most journals and scientists consider an alpha larger than 0.05 untrustworthy. We therefore recommend using alpha values smaller than or queal to 0.05'
            )
        }

        # test wether switching have been analyzed
        if (!any(!is.na(
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value
        ))) {
            stop(
                'The analsis of isoform switching must be performed before functional consequences can be analyzed. Please run ?detectIsoformSwitching and try again.'
            )
        }

        acceptedTypes <- c(
            # Transcript
            'tss',
            'tts',
            'last_exon',
            'isoform_length',
            'exon_number',
            'intron_structure',
            'intron_retention',
            'isoform_class_code',
            # cpat
            'coding_potential',
            # ORF
            'ORF_genomic',
            'ORF_length',
            '5_utr_length',
            '3_utr_length',
            # seq similarity
            'isoform_seq_similarity',
            'ORF_seq_similarity',
            '5_utr_seq_similarity',
            '3_utr_seq_similarity',
            # ORF
            'NMD_status',
            # pfam
            'domains_identified',
            'genomic_domain_position',
            'domain_length',
            # SignalIP
            'signal_peptide_identified',
            # IDR
            'IDR_identified',
            'IDR_length',
            'IDR_type',
            # DeepLoc3
            'sub_cell_location',
            'solubility_status'
        )

        if (!all(consequencesToAnalyze %in% c('all', acceptedTypes))) {
            stop(
                paste(
                    'The argument(s) supplied to \'typeOfconsequence\' are not accepted.',
                    'Please see ?analyzeSwitchConsequences for description of which strings are allowed.',
                    'The problem is:',
                    paste(setdiff(
                        consequencesToAnalyze , c('all', acceptedTypes)
                    ), collapse = ', '),
                    sep = ' '
                )
            )
        }

        if ('all' %in% consequencesToAnalyze) {
            consequencesToAnalyze <- acceptedTypes
        }

        ## Test whether annotation is advailable
        if ('intron_retention'  %in% consequencesToAnalyze) {
            if (is.null(switchAnalyzeRlist$intronRetentionAnalysis) & is.null( switchAnalyzeRlist$AlternativeSplicingAnalysis)) {
                stop(
                    'To test for intron retention alternative splicing must first be classified. Please run analyzeIntronRetention() and try again.'
                )
            }
        }

        if (grepl('cufflinks', switchAnalyzeRlist$sourceId)) {
            if ('isoform_class_code'  %in% consequencesToAnalyze) {
                if (!'class_code' %in%
                    colnames(switchAnalyzeRlist$isoformFeatures)
                ) {
                    stop(
                        'The switchAnalyzeRlist does not contail the calss_code information'
                    )
                }
            }
        } else {
            consequencesToAnalyze <-
                setdiff(consequencesToAnalyze, 'isoform_class_code')
        }

        if (any(c('ORF_genomic', 'ORF_length', 'NMD_status') %in%
                consequencesToAnalyze
        )) {
            if (!'PTC' %in% colnames(switchAnalyzeRlist$isoformFeatures)) {
                stop(
                    'To test differences in ORFs or PCT, ORF must be annotated. Please run annotateORF() and try again'
                )
            }
        }
        if ('coding_potential'  %in% consequencesToAnalyze) {
            if (!'codingPotential' %in%
                colnames(switchAnalyzeRlist$isoformFeatures)
            ) {
                stop(
                    'To test differences in coding_potential, the result of the CPAT analysis must be advailable. Please run addCPATanalysis() and try again'
                )
            }
        }
        if (any(
            c(
                'domains_identified',
                'domain_length',
                'genomic_domain_position'
            )  %in% consequencesToAnalyze
        )) {
            if (is.null(switchAnalyzeRlist$domainAnalysis)) {
                stop(
                    'To test differences in protein domains, the result of the Pfam analysis must be advailable. Please run addPFAManalysis() and try again'
                )
            }
        }
        if ('signal_peptide_identified'  %in% consequencesToAnalyze) {
            if (is.null(switchAnalyzeRlist$signalPeptideAnalysis)) {
                stop(
                    'To test differences in signal peptides, the result of the SignalP analysis must be advailable. Please run addSignalIPanalysis() and try again'
                )
            }
        }

        if ('sub_cell_location'  %in% consequencesToAnalyze) {
            if ( ! 'sub_cell_location' %in% colnames(switchAnalyzeRlist$isoformFeatures)) {
                stop(
                    'Cannot test for differences in sub-cellular location as such results are not annotated'
                )
            }
        }
        if ('solubility_status'  %in% consequencesToAnalyze) {
            if ( ! 'solubility_status' %in% colnames(switchAnalyzeRlist$isoformFeatures)) {
                stop(
                    'Cannot test for differences in solubility status as such results are not annotated'
                )
            }
        }


        if (!is.numeric(ntCutoff)) {
            stop('The ntCutoff arugment must be an numeric')
        }
        if (ntCutoff <= 0) {
            stop('The ntCutoff arugment must be an numeric > 0')
        }

        if (!is.null(ntFracCutoff)) {
            if (ntFracCutoff <= 0 | ntFracCutoff > 1) {
                stop(
                    'The ntFracCutoff arugment must be a numeric in the interval (0,1]. Use NULL to disable.'
                )
            }
        }

        ### test sequence annotation
        if (any(
            consequencesToAnalyze %in% c(
                'isoform_seq_similarity',
                '5_utr_seq_similarity',
                '3_utr_seq_similarity'
            )
        )) {
            if (!any(names(switchAnalyzeRlist) == 'ntSequence')) {
                stop(
                    'The transcrip nucleotide sequences must be added to the switchAnalyzeRlist before overlap analysis can be performed. Please run \'extractSequence\' and try again.'
                )
            }
        }
        if (any(consequencesToAnalyze %in% c('ORF_seq_similarity'))) {
            if (!any(names(switchAnalyzeRlist) == 'aaSequence')) {
                stop(
                    'The transcrip ORF amino acid sequences must be added to the switchAnalyzeRlist before ORF overlap analysis can be performed. Please run \'extractSequence\' and try again.'
                )
            }
        }

        if( 'sub_cell_location' %in% consequencesToAnalyze) {
            if( ! 'sub_cell_location' %in% colnames(switchAnalyzeRlist$isoformFeatures)) {
                'Subcellular localizations have not been analyzed. Please annotate and try again.'
            }
        }
        if( 'solubility_status' %in% consequencesToAnalyze) {
            if( ! 'solubility_status' %in% colnames(switchAnalyzeRlist$isoformFeatures)) {
                'Subcellular localizations have not been analyzed. Please annotate and try again.'
            }
        }

    }

    ### Extract and massage data
    if (TRUE) {
        if (!is.null(ntFracCutoff)) {
            fractionFilter <- TRUE
        } else {
            fractionFilter <- FALSE
        }

        isoformsToAnalyze <-  c(upIso, downIso)
        names(isoformsToAnalyze) <- c('up', 'down')

        ### transcript structure
        exonData <-
            switchAnalyzeRlist$exons[which(
                switchAnalyzeRlist$exons$isoform_id %in% isoformsToAnalyze
                ),
                'isoform_id'
            ] # tss, tts, isoform_length, exon_number
        exonDataList <- split(exonData, f = exonData$isoform_id)

        ### transcript annotation
        columnsToExtract <- 'isoform_id'
        if ('isoform_class_code' %in% consequencesToAnalyze) {
            columnsToExtract <- c(columnsToExtract, 'class_code')
        }
        if ('intron_retention'  %in% consequencesToAnalyze) {
            columnsToExtract <- c(columnsToExtract, 'IR')
        }
        if ('NMD_status'         %in% consequencesToAnalyze) {
            columnsToExtract <- c(columnsToExtract, 'PTC')
        }
        if ('coding_potential'   %in% consequencesToAnalyze) {
            columnsToExtract <- c(columnsToExtract, 'codingPotential')
        }
        if ('sub_cell_location'   %in% consequencesToAnalyze) {
            columnsToExtract <- c(columnsToExtract, 'sub_cell_location')
        }
        if ('solubility_status'   %in% consequencesToAnalyze) {
            columnsToExtract <- c(columnsToExtract, 'solubility_status')
        }

        transcriptData <-
            unique(switchAnalyzeRlist$isoformFeatures[
                which(
                    switchAnalyzeRlist$isoformFeatures$isoform_id %in%
                        isoformsToAnalyze
                ), columnsToExtract,
                drop =FALSE
            ])

        ### ORF data
        columnsToExtract2 <- 'isoform_id'
        if ('ORF_genomic'  %in% consequencesToAnalyze) {
            columnsToExtract2 <-
                c(columnsToExtract2,
                  c('orfStartGenomic', 'orfEndGenomic'))
        }
        if (any(
            c(
                'ORF_length',
                '5_utr_length',
                '3_utr_length',
                'ORF_seq_similarity',
                '5_utr_seq_similarity',
                '3_utr_seq_similarity',
                'domains_identified',
                'genomic_domain_position',
                'domain_length',
                'signal_peptide_identified',
                'IDR_identified',
                'IDR_length',
                'IDR_type',
                'solubility_status',
                'sub_cell_location'
            ) %in% consequencesToAnalyze
        )) {
            columnsToExtract2 <-
                c(
                    columnsToExtract2,
                    c(
                        'orfTransciptStart',
                        'orfTransciptEnd',
                        'orfTransciptLength'
                    )
                )
        }
        if (length(columnsToExtract2) > 1) {
            orfData <- unique(switchAnalyzeRlist$orfAnalysis[
                which(
                    switchAnalyzeRlist$orfAnalysis$isoform_id %in%
                        isoformsToAnalyze
                    ),
                columnsToExtract2
            ])

            transcriptData <-
                dplyr::left_join(transcriptData,
                      orfData,
                      by = 'isoform_id',
                      all.x = TRUE)
        }

        ### intron retention data
        if ('intron_retention' %in% consequencesToAnalyze) {
            if( is.null( switchAnalyzeRlist$AlternativeSplicingAnalysis) ) {
                localIRdata <-
                    switchAnalyzeRlist$intronRetentionAnalysis[which(
                        switchAnalyzeRlist$intronRetentionAnalysis$isoform_id %in%
                            isoformsToAnalyze
                    ), ]
            } else {
                localIRdata <-
                    switchAnalyzeRlist$AlternativeSplicingAnalysis[which(
                        switchAnalyzeRlist$AlternativeSplicingAnalysis$isoform_id %in%
                            isoformsToAnalyze
                    ), ]
            }

            if (nrow(localIRdata) != 2) {
                warning(
                    paste(
                        'There was a problem with the extraction if intron retentions -',
                        'please contact the developers with this example so they can fix it.',
                        'For now they are ignored. The isoforms affected are:',
                        paste(isoformsToAnalyze, collapse = ', '),
                        sep = ' '
                    )
                )
                consequencesToAnalyze <-
                    consequencesToAnalyze[which(
                        !consequencesToAnalyze %in% c('intron_retention')
                    )]
            } else {
                localIRdata$irCoordinats <-
                    paste(localIRdata$IR_genomic_start,
                          localIRdata$IR_genomic_end,
                          sep = ':')
            }
        }

        ### If nessesary extract data to remove
        onPlusStrand <- as.character(strand(exonData)[1]) == '+'
        if( 'wasTrimmed' %in% colnames(switchAnalyzeRlist$orfAnalysis) ) {

            if(onPlusStrand) {
                localTrimmed <- switchAnalyzeRlist$orfAnalysis[
                    which(switchAnalyzeRlist$orfAnalysis$isoform_id %in% isoformsToAnalyze),
                    c('isoform_id','wasTrimmed','trimmedStartGenomic','orfEndGenomic')
                ]
            } else {
                localTrimmed <- switchAnalyzeRlist$orfAnalysis[
                    which(switchAnalyzeRlist$orfAnalysis$isoform_id %in% isoformsToAnalyze),
                    c('isoform_id','wasTrimmed','trimmedStartGenomic','orfStartGenomic')
                    ]
                colnames(localTrimmed) <- c('isoform_id','wasTrimmed','trimmedStartGenomic','orfEndGenomic')
            }

            if(any(localTrimmed$wasTrimmed, na.rm = TRUE)) {
                localTrimmed <- localTrimmed[which(
                    localTrimmed$wasTrimmed
                ),]

                regionToOmmit <- GenomicRanges::reduce(IRanges(localTrimmed$trimmedStartGenomic, localTrimmed$orfEndGenomic))
            }
        }

        ### domain data
        if (any(
            c(
                'domains_identified',
                'genomic_domain_position',
                'domain_length'
            )  %in% consequencesToAnalyze
        )) {
            domanData <-
                switchAnalyzeRlist$domainAnalysis[which(
                    switchAnalyzeRlist$domainAnalysis$isoform_id %in%
                        isoformsToAnalyze
                ), ]
            domanData$isoform_id <-
                factor(domanData$isoform_id, levels = isoformsToAnalyze)
            domanDataSplit <-
                split(domanData[, c(
                    'hmm_name',
                    'pfamStartGenomic',
                    'pfamEndGenomic',
                    'orf_aa_start',
                    'orf_aa_end'
                )], f = domanData$isoform_id)

            ### Remove those overlapping trimmed regions
            if( exists('regionToOmmit') ) {
                domanDataSplit <- lapply(domanDataSplit, function(aSet) { # aSet <- domanDataSplit[[2]]
                    if( onPlusStrand ) {
                        aSet[which(
                            ! overlapsAny(
                                IRanges(aSet$pfamStartGenomic, aSet$pfamEndGenomic),
                                regionToOmmit
                            )
                        ),]
                    } else {
                        aSet[which(
                            ! overlapsAny(
                                IRanges(aSet$pfamEndGenomic, aSet$pfamStartGenomic),
                                regionToOmmit
                            )
                        ),]
                    }

                })
            }

            # overwrite if no ORF is detected
            isNAnames <-
                transcriptData$isoform_id[which(
                    is.na(transcriptData$orfTransciptLength)
                )]
            if (length(isNAnames)) {
                isNAindex <- which(names(domanDataSplit) %in% isNAnames)
                domanDataSplit[isNAindex] <-
                    lapply(domanDataSplit[isNAindex], function(aDF) {
                        aDF[0, ]
                    })
            }
        }

        if ( any( c('IDR_identified','IDR_type','IDR_length') %in% consequencesToAnalyze ) ) {
            idrData <-
                switchAnalyzeRlist$idrAnalysis[which(
                    switchAnalyzeRlist$idrAnalysis$isoform_id %in%
                        isoformsToAnalyze
                ), ]
            idrData$isoform_id <-
                factor(idrData$isoform_id, levels = isoformsToAnalyze)
            idrDataSplit <-
                split(idrData[, c(
                    'idrStartGenomic',
                    'idrEndGenomic',
                    'orf_aa_start',
                    'orf_aa_end',
                    'idr_type'
                )], f = idrData$isoform_id)

            ### Remove those overlapping trimmed regions
            if( exists('regionToOmmit') ) {
                idrDataSplit <- lapply(idrDataSplit, function(aSet) { # aSet <- idrDataSplit[[2]]
                    if( onPlusStrand ) {
                        aSet[which(
                            ! overlapsAny(
                                IRanges(aSet$idrStartGenomic, aSet$idrEndGenomic),
                                regionToOmmit
                            )
                        ),]
                    } else {
                        aSet[which(
                            ! overlapsAny(
                                IRanges(aSet$idrEndGenomic, aSet$idrStartGenomic),
                                regionToOmmit
                            )
                        ),]
                    }
                })
            }

            ### overwrite if no ORF is detected
            isNAnames <-
                transcriptData$isoform_id[which(
                    is.na(transcriptData$orfTransciptLength)
                )]
            if (length(isNAnames)) {
                isNAindex <- which(names(idrDataSplit) %in% isNAnames)
                idrDataSplit[isNAindex] <-
                    lapply(idrDataSplit[isNAindex], function(aDF) {
                        aDF[0, ]
                    })
            }
        }

        ### peptide data
        if ('signal_peptide_identified'  %in% consequencesToAnalyze) {
            peptideData <-
                switchAnalyzeRlist$signalPeptideAnalysis[which(
                    switchAnalyzeRlist$signalPeptideAnalysis$isoform_id %in%
                        isoformsToAnalyze
                ), ]
            peptideData$isoform_id <-
                factor(peptideData$isoform_id, levels = isoformsToAnalyze)
            peptideDataSplit <-
                split(peptideData, f = peptideData$isoform_id)


            ### Remove those overlapping trimmed regions
            if( exists('regionToOmmit') ) {
                peptideDataSplit <- lapply(peptideDataSplit, function(aSet) { # aSet <- peptideDataSplit[[2]]
                    aSet[which(
                        ! overlapsAny(
                            IRanges(aSet$genomicClevageAfter, aSet$genomicClevageAfter),
                            regionToOmmit
                        )
                    ),]
                })
            }


            # overwrite if no ORF is detected
            isNAnames <-
                transcriptData$isoform_id[which(
                    is.na(transcriptData$orfTransciptLength)
                )]
            if (length(isNAnames)) {
                isNAindex <- which(names(peptideDataSplit) %in% isNAnames)
                peptideDataSplit[isNAindex] <-
                    lapply(peptideDataSplit[isNAindex], function(aDF) {
                        aDF[0, ]
                    })
            }
        }

        ### Extract isoform length
        if (any(
            c(
                'isoform_length',
                '5_utr_length',
                '3_utr_length',
                'ORF_length',
                'isoform_seq_similarity',
                '5_utr_seq_similarity',
                '3_utr_seq_similarity',
                'ORF_seq_similarity'
            ) %in% consequencesToAnalyze
        )) {
            isoform_length <- sapply(exonDataList, function(x)
                sum(width(x)))
            transcriptData$length <-
                isoform_length[match(transcriptData$isoform_id,
                                     names(isoform_length))]
        }

        ### Sequence overlap
        if (any(
            c(
                'isoform_seq_similarity',
                '5_utr_seq_similarity',
                '3_utr_seq_similarity'
            ) %in% consequencesToAnalyze
        )) {
            ntSeq <-
                switchAnalyzeRlist$ntSequence[which(
                    names(switchAnalyzeRlist$ntSequence) %in% c(upIso, downIso)
                )]
            if (length(ntSeq) == 2) {
                upNtSeq   <- ntSeq[upIso]
                downNtSeq <- ntSeq[downIso]
            } else {
                #warning('The amino acid sequence of some of the transcripts to compare is not pressent - \'isoform_seq_similarity\',\'5_utr_seq_similarity\' and \'3_utr_seq_similarity\' will therefor be ignored.')
                consequencesToAnalyze <-
                    consequencesToAnalyze[which(
                        !consequencesToAnalyze %in% c(
                            'isoform_seq_similarity',
                            '5_utr_seq_similarity',
                            '3_utr_seq_similarity'
                        )
                    )]
            }
        }
        if ('ORF_seq_similarity' %in% consequencesToAnalyze) {
            aaSeq <-
                switchAnalyzeRlist$aaSequence[which(
                    names(switchAnalyzeRlist$aaSequence) %in%
                        transcriptData$isoform_id[which(
                            !is.na(transcriptData$orfTransciptLength)
                        )]
                )]
            if (length(aaSeq) == 2) {
                upAAseq   <- aaSeq[upIso]
                downAAseq <- aaSeq[downIso]
            }
        }

        if (length(consequencesToAnalyze) == 0) {
            return(NULL)
        }
        if (nrow(transcriptData) != 2) {
            return(NULL)
        }

    }

    ### Analyze differences
    if (TRUE) {
        ### Make result data.frame
        isoComparison <-
            data.frame(
                isoformUpregulated = upIso,
                isoformDownregulated = downIso,
                featureCompared = consequencesToAnalyze,
                isoformsDifferent = NA
            )
        if (addDescription) {
            isoComparison$switchConsequence <- NA
        }

        ### analyze  isoform differnces
        if ('tss'                       %in% consequencesToAnalyze) {
            localExonData <- unlist(range(exonDataList))

            if (as.character(localExonData@strand[1]) == '+') {
                tssCoordinats <- start(localExonData)

                ### test
                tssDifferent <-
                    abs(tssCoordinats[1] - tssCoordinats[2]) > ntCutoff
                mostUpstream <-
                    names(localExonData)[which.min(tssCoordinats)]
            } else {
                tssCoordinats <- end(localExonData)

                # test
                tssDifferent <-
                    abs(tssCoordinats[1] - tssCoordinats[2]) > ntCutoff
                mostUpstream <-
                    names(localExonData)[which.max(tssCoordinats)]
            }

            # make repport
            localIndex <-
                which(isoComparison$featureCompared == 'tss')
            isoComparison$isoformsDifferent[localIndex] <-
                tssDifferent

            # add description
            if (tssDifferent & addDescription) {
                switchMoreUpstram <- mostUpstream == upIso

                if (switchMoreUpstram) {
                    isoComparison$switchConsequence[localIndex] <-
                        'Tss more upstream'
                } else {
                    isoComparison$switchConsequence[localIndex] <-
                        'Tss more downstream'
                }
            }
        }

        if ('tts'                       %in% consequencesToAnalyze) {
            localExonData <- unlist(range(exonDataList))

            if (as.character(localExonData@strand[1]) == '+') {
                ttsCoordinats <- end(localExonData)

                # test
                ttsDifferent <-
                    abs(ttsCoordinats[1] - ttsCoordinats[2]) > ntCutoff
                mostDownstream <-
                    names(localExonData)[which.max(ttsDifferent)]
            } else {
                ttsCoordinats <- start(localExonData)

                # test
                ttsDifferent <-
                    abs(ttsCoordinats[1] - ttsCoordinats[2]) > ntCutoff
                mostDownstream <-
                    names(localExonData)[which.min(ttsDifferent)]
            }

            # make repport
            localIndex <-
                which(isoComparison$featureCompared == 'tts')
            isoComparison$isoformsDifferent[localIndex] <-
                ttsDifferent

            # add description
            if (ttsDifferent & addDescription) {
                switchMoreDownstream <- mostDownstream == upIso

                if (switchMoreDownstream) {
                    isoComparison$switchConsequence[localIndex] <-
                        'Tts more downstream'
                } else {
                    isoComparison$switchConsequence[localIndex] <-
                        'Tts more upstream'
                }
            }
        }

        if ('last_exon'                 %in% consequencesToAnalyze) {
            isPlusStrand <- as.logical(exonData@strand[1] == '+')
            if (isPlusStrand) {
                lastExons <- lapply(exonDataList, function(x)
                    x[length(x), 0])
            } else {
                lastExons <- lapply(exonDataList, function(x)
                    x[1        , 0])
            }

            lastExonDifferent <-
                !overlapsAny(lastExons[[1]], lastExons[[2]])

            # make repport
            localIndex <-
                which(isoComparison$featureCompared == 'last_exon')
            isoComparison$isoformsDifferent[localIndex] <-
                lastExonDifferent

            # add description
            if (lastExonDifferent & addDescription) {
                if (isPlusStrand) {
                    localEndCoordinats <-
                        sapply(lastExons, function(aGRange)
                            end(aGRange))
                    mostDownstream <-
                        names(localEndCoordinats)[which.max(localEndCoordinats)]
                } else {
                    localEndCoordinats <-
                        sapply(lastExons, function(aGRange)
                            start(aGRange))
                    mostDownstream <-
                        names(localEndCoordinats)[which.min(localEndCoordinats)]
                }

                switchMoreDownstream <- mostDownstream == upIso

                if (switchMoreDownstream) {
                    isoComparison$switchConsequence[localIndex] <-
                        'Last exon more downstream'
                } else {
                    isoComparison$switchConsequence[localIndex] <-
                        'Last exon more upstream'
                }
            }
        }

        if ('isoform_length'            %in% consequencesToAnalyze) {
            differentLength <- abs(diff(isoform_length)) > ntCutoff

            if (fractionFilter) {
                downLength <-
                    isoform_length[which(names(isoform_length) == downIso)]
                fractionDifference <-
                    abs(diff(isoform_length)) / downLength > ntFracCutoff

                differentLength <-
                    differentLength & fractionDifference
            }

            # make repport
            localIndex <-
                which(isoComparison$featureCompared == 'isoform_length')
            isoComparison$isoformsDifferent[localIndex] <-
                differentLength

            # make description
            if (differentLength & addDescription) {
                lengthGain <-
                    names(isoform_length)[which.max(isoform_length)] == upIso

                if (lengthGain) {
                    isoComparison$switchConsequence[localIndex] <- 'Length gain'
                } else {
                    isoComparison$switchConsequence[localIndex] <- 'Length loss'
                }

            }
        }

        if ('isoform_seq_similarity'    %in% consequencesToAnalyze) {
            localAlignment <-
                Biostrings::pairwiseAlignment(pattern = upNtSeq,
                                  subject = downNtSeq,
                                  type = 'global')

            overlapSize <- min(c(nchar(gsub(
                '-', '', as.character(localAlignment@subject)
            )),
            nchar(gsub(
                '-', '', as.character(localAlignment@pattern)
            ))))
            totalWidth <- width(localAlignment@subject@unaligned) +
                width(localAlignment@pattern@unaligned) -
                overlapSize

            jcDist <- overlapSize / totalWidth

            differentIsoformOverlap <-
                jcDist < ntJCsimCutoff & totalWidth - overlapSize > ntCutoff

            # make repport
            localIndex <-
                which(isoComparison$featureCompared == 'isoform_seq_similarity')
            isoComparison$isoformsDifferent[localIndex] <-
                differentIsoformOverlap

            # make description
            if (differentIsoformOverlap & addDescription) {
                lengthGain <-
                    names(isoform_length)[which.max(isoform_length)] == upIso

                if (lengthGain) {
                    isoComparison$switchConsequence[localIndex] <- 'Length gain'
                } else {
                    isoComparison$switchConsequence[localIndex] <- 'Length loss'
                }

            }
        }

        if ('exon_number'               %in% consequencesToAnalyze) {
            localNrExons <- sapply(exonDataList, length)
            exonDifferent <- localNrExons[1] != localNrExons[2]

            # make repport
            localIndex <-
                which(isoComparison$featureCompared == 'exon_number')
            isoComparison$isoformsDifferent[localIndex] <-
                exonDifferent

            # make description
            if (exonDifferent & addDescription) {
                switchGainsExons <-
                    names(localNrExons)[which.max(localNrExons)] == upIso

                if (switchGainsExons) {
                    isoComparison$switchConsequence[localIndex] <- 'Exon gain'
                } else {
                    isoComparison$switchConsequence[localIndex] <- 'Exon loss'
                }
            }

        }

        if ('intron_structure'          %in% consequencesToAnalyze) {
            localIntrons <-
                lapply(exonDataList, function(aGR) {
                    gaps(ranges(aGR))
                })

            differentintron_structure <- !any(
                all(localIntrons[[1]] %in% localIntrons[[2]]),
                all(localIntrons[[2]] %in% localIntrons[[1]])
            )

            # make repport
            localIndex <-
                which(isoComparison$featureCompared == 'intron_structure')
            isoComparison$isoformsDifferent[localIndex] <-
                differentintron_structure
        }

        if ('intron_retention'          %in% consequencesToAnalyze) {
            if (all(!is.na(transcriptData$IR))) {
                differentNrIR <-
                    localIRdata$irCoordinats[1] != localIRdata$irCoordinats[2]

                # Make repport
                localIndex <-
                    which(isoComparison$featureCompared == 'intron_retention')
                isoComparison$isoformsDifferent[localIndex] <-
                    differentNrIR

                # make description
                if (differentNrIR & addDescription) {
                    if (transcriptData$IR[1] == transcriptData$IR[2]) {
                        isoComparison$switchConsequence[localIndex] <-
                            'Intron retention switch'
                    } else if (transcriptData$isoform_id[which.max(
                        transcriptData$IR
                    )] == upIso) {
                        isoComparison$switchConsequence[localIndex] <-
                            'Intron retention gain'
                    } else {
                        isoComparison$switchConsequence[localIndex] <-
                            'Intron retention loss'
                    }
                }
            }
        }

        if ('isoform_class_code'        %in% consequencesToAnalyze) {
            differentAnnotation <-
                transcriptData$class_code[1] != transcriptData$class_code[2]

            # make repport
            localIndex <-
                which(isoComparison$featureCompared == 'isoform_class_code')
            isoComparison$isoformsDifferent[localIndex] <-
                differentAnnotation
        }

        if ('ORF_genomic'               %in% consequencesToAnalyze) {
            localIndex <- which(isoComparison$featureCompared == 'ORF_genomic')

            # if both have ORFs
            if (all(!is.na(transcriptData$orfStartGenomic))) {
                genomicOrfDifferent <-
                    transcriptData$orfStartGenomic[1] !=
                    transcriptData$orfStartGenomic[2] |
                    transcriptData$orfEndGenomic[1] !=
                    transcriptData$orfEndGenomic[2]

                # make repport
                isoComparison$isoformsDifferent[localIndex] <-
                    genomicOrfDifferent
                # If only 1 ORF
            } else if (sum(!is.na(transcriptData$orfStartGenomic)) == 1) {
                # make repport
                isoComparison$isoformsDifferent[localIndex] <- TRUE

            } else {
                isoComparison$isoformsDifferent[localIndex] <- FALSE
            }
        }

        if ('ORF_length'                %in% consequencesToAnalyze) {
            localIndex <- which(isoComparison$featureCompared == 'ORF_length')
            # if both ORFs are annotated
            if (all(!is.na(transcriptData$orfTransciptLength))) {
                orfDifferent <-
                    abs(diff(transcriptData$orfTransciptLength)) > ntCutoff

                if (fractionFilter) {
                    downLength <-
                        isoform_length[which(names(isoform_length) == downIso)]
                    fractionDifference <-
                        abs(diff(transcriptData$orfTransciptLength)) /
                        downLength > ntFracCutoff

                    orfDifferent <-
                        orfDifferent & fractionDifference
                }

                # make repport
                isoComparison$isoformsDifferent[localIndex] <-
                    orfDifferent

                # make description
                if (orfDifferent & addDescription) {
                    orfGain <-
                        transcriptData$isoform_id[which.max(
                            transcriptData$orfTransciptLength
                        )] == upIso

                    if (orfGain) {
                        isoComparison$switchConsequence[localIndex] <-
                            'ORF is longer'
                    } else {
                        isoComparison$switchConsequence[localIndex] <-
                            'ORF is shorter'
                    }
                }
                #If only one ORF is annotated
            } else if (sum(!is.na(transcriptData$orfTransciptLength)) == 1) {
                # make repport
                isoComparison$isoformsDifferent[localIndex] <- TRUE

                if (addDescription) {
                    orfLoss <-
                        is.na(transcriptData$orfTransciptLength[which(
                            transcriptData$isoform_id == upIso
                        )])

                    if (orfLoss) {
                        isoComparison$switchConsequence[localIndex] <-
                            'Complete ORF loss'
                    } else {
                        isoComparison$switchConsequence[localIndex] <-
                            'Complete ORF gain'
                    }
                }
            } else {
                isoComparison$isoformsDifferent[localIndex] <- FALSE
            }
        }

        if ('ORF_seq_similarity'        %in% consequencesToAnalyze) {
            localIndex <-
                which(isoComparison$featureCompared == 'ORF_seq_similarity')
            # if both have ORFs
            if (all(!is.na(transcriptData$orfTransciptLength))) {
                if (length(aaSeq) == 2) {
                    localAlignment <-
                        Biostrings::pairwiseAlignment(
                            pattern = upAAseq,
                            subject = downAAseq,
                            type = 'global'
                        )

                    overlapSize <- min(c(nchar(
                        gsub(
                            '-',
                            '',
                            as.character(localAlignment@subject)
                        )
                    ),
                    nchar(
                        gsub(
                            '-',
                            '',
                            as.character(localAlignment@pattern)
                        )
                    )))
                    totalWidth <-
                        width(localAlignment@subject@unaligned) +
                        width(localAlignment@pattern@unaligned) -
                        overlapSize

                    jcDist <- overlapSize / totalWidth

                    differentORFoverlap <-
                        jcDist < AaJCsimCutoff &
                        (totalWidth - overlapSize) > AaCutoff

                    # make repport
                    isoComparison$isoformsDifferent[localIndex] <-
                        differentORFoverlap

                    # make description
                    if (differentORFoverlap & addDescription) {
                        lengthGain <-
                            transcriptData$isoform_id[which.max(
                                transcriptData$orfTransciptLength
                            )] == upIso

                        if (lengthGain) {
                            isoComparison$switchConsequence[localIndex] <-
                                'ORF is longer'
                        } else {
                            isoComparison$switchConsequence[localIndex] <-
                                'ORF is shorter'
                        }

                    }
                }
                # If only 1 ORF is annotated
            } else if (sum(!is.na(transcriptData$orfTransciptLength)) == 1) {
                # make repport
                isoComparison$isoformsDifferent[localIndex] <- TRUE

                if (addDescription) {
                    orfLoss <-
                        is.na(transcriptData$orfTransciptLength[which(
                            transcriptData$isoform_id == upIso
                        )])

                    if (orfLoss) {
                        isoComparison$switchConsequence[localIndex] <-
                            'Complete ORF loss'
                    } else {
                        isoComparison$switchConsequence[localIndex] <-
                            'Complete ORF gain'
                    }
                }
            } else {
                isoComparison$isoformsDifferent[localIndex] <- FALSE
            }
        }

        if ('5_utr_length'              %in% consequencesToAnalyze) {
            if (all(!is.na(transcriptData$orfTransciptStart))) {
                utr5Different <-
                    abs(diff(transcriptData$orfTransciptStart)) > ntCutoff

                if (fractionFilter) {
                    downLength <-
                        isoform_length[which(names(isoform_length) == downIso)]
                    fractionDifference <-
                        abs(diff(transcriptData$orfTransciptStart)) /
                        downLength > ntFracCutoff

                    utr5Different <-
                        utr5Different & fractionDifference
                }


                # make repport
                localIndex <-
                    which(isoComparison$featureCompared == '5_utr_length')
                isoComparison$isoformsDifferent[localIndex] <-
                    utr5Different

                # make description
                if (utr5Different & addDescription) {
                    utr5Gain <-
                        transcriptData$isoform_id[which.max(
                            transcriptData$orfTransciptStart
                        )] == upIso

                    if (utr5Gain) {
                        isoComparison$switchConsequence[localIndex] <-
                            '5UTR is longer'
                    } else {
                        isoComparison$switchConsequence[localIndex] <-
                            '5UTR is shorter'
                    }
                }
            }
        }

        if ('5_utr_seq_similarity'      %in% consequencesToAnalyze) {
            if (all(!is.na(transcriptData$orfTransciptStart))) {
                localUpNt   <-
                    subseq(upNtSeq,
                           1,
                           transcriptData$orfTransciptStart[which(
                               transcriptData$isoform_id == upIso
                           )] - 1)
                localDownNt <-
                    subseq(downNtSeq,
                           1,
                           transcriptData$orfTransciptStart[which(
                               transcriptData$isoform_id == downIso
                           )] -1)

                # if one of them have no UTR
                if (width(localUpNt) > 0 &
                    width(localDownNt) > 0) {
                    localAlignment <-
                        Biostrings::pairwiseAlignment(
                            pattern = localUpNt,
                            subject = localDownNt,
                            type = 'overlap'
                        )

                    overlapSize <- min(c(nchar(
                        gsub(
                            '-',
                            '',
                            as.character(localAlignment@subject)
                        )
                    ),
                    nchar(
                        gsub(
                            '-',
                            '',
                            as.character(localAlignment@pattern)
                        )
                    )))
                    totalWidth <-
                        width(localAlignment@subject@unaligned) +
                        width(localAlignment@pattern@unaligned) -
                        overlapSize

                    jcDist <- overlapSize / totalWidth
                } else {
                    overlapSize <- 0
                    totalWidth <-
                        abs(width(localUpNt) - width(localDownNt))
                    jcDist <- 0
                }

                differenttUTRoverlap <-
                    jcDist < ntJCsimCutoff & totalWidth - overlapSize > ntCutoff

                # make repport
                localIndex <-
                    which(
                        isoComparison$featureCompared == '5_utr_seq_similarity'
                    )
                isoComparison$isoformsDifferent[localIndex] <-
                    differenttUTRoverlap

                # make description
                if (differenttUTRoverlap & addDescription) {
                    utr5Gain <-
                        transcriptData$isoform_id[which.max(
                            transcriptData$orfTransciptStart
                        )] == upIso

                    if (utr5Gain) {
                        isoComparison$switchConsequence[localIndex] <-
                            '5UTR is longer'
                    } else {
                        isoComparison$switchConsequence[localIndex] <-
                            '5UTR is shorter'
                    }

                }
            }
        }

        if ('3_utr_length'              %in% consequencesToAnalyze) {
            if (all(!is.na(transcriptData$orfTransciptEnd))) {
                transcriptData$utr3length <-
                    transcriptData$length - (transcriptData$orfTransciptEnd + 1)

                utr3Different <-
                    abs(diff(transcriptData$utr3length)) > ntCutoff

                if (fractionFilter) {
                    downLength <-
                        isoform_length[which(names(isoform_length) == downIso)]
                    fractionDifference <-
                        abs(diff(transcriptData$utr3length)) /
                        downLength > ntFracCutoff

                    utr3Different <-
                        utr3Different & fractionDifference
                }

                # make repport
                localIndex <-
                    which(isoComparison$featureCompared == '3_utr_length')
                isoComparison$isoformsDifferent[localIndex] <-
                    utr3Different

                # make description
                if (utr3Different & addDescription) {
                    utr3Gain <-
                        transcriptData$isoform_id[which.max(
                            transcriptData$utr3length
                        )] == upIso

                    if (utr3Gain) {
                        isoComparison$switchConsequence[localIndex] <-
                            '3UTR is longer'
                    } else {
                        isoComparison$switchConsequence[localIndex] <-
                            '3UTR is shorter'
                    }
                }
            }
        }

        if ('3_utr_seq_similarity'      %in% consequencesToAnalyze) {
            if (all(!is.na(transcriptData$orfTransciptEnd))) {
                transcriptData$utr3length <-
                    transcriptData$length - (transcriptData$orfTransciptEnd + 1)

                up3UTRstart   <-
                    transcriptData$orfTransciptEnd[which(
                        transcriptData$isoform_id == upIso
                    )] + 4 # +4 because the stop codon is not included in the ORF
                down3UTRstart <-
                    transcriptData$orfTransciptEnd[which(
                        transcriptData$isoform_id == downIso
                    )] + 4 # +4 because the stop codon is not included in the ORF

                # if 3UTR somehow is annotated as longer set it to the length +1 (then the sequence will be length 0)
                if (up3UTRstart   > width(upNtSeq)   + 1) {
                    up3UTRstart   <- width(upNtSeq)  + 1
                }
                if (down3UTRstart > width(downNtSeq) + 1) {
                    down3UTRstart <- width(downNtSeq) + 1
                }

                localUpNt   <-
                    subseq(upNtSeq,   up3UTRstart,   width(upNtSeq))
                localDownNt <-
                    subseq(downNtSeq, down3UTRstart, width(downNtSeq))

                ### if one of them have no UTR
                if (width(localUpNt) > 0 &
                    width(localDownNt) > 0) {
                    localAlignment <-
                        Biostrings::pairwiseAlignment(
                            pattern = localUpNt,
                            subject = localDownNt,
                            type = 'overlap'
                        )

                    overlapSize <- min(c(nchar(
                        gsub(
                            '-',
                            '',
                            as.character(localAlignment@subject)
                        )
                    ),
                    nchar(
                        gsub(
                            '-',
                            '',
                            as.character(localAlignment@pattern)
                        )
                    )))
                    totalWidth <-
                        width(localAlignment@subject@unaligned) +
                        width(localAlignment@pattern@unaligned) -
                        overlapSize

                    jcDist <- overlapSize / totalWidth
                } else {
                    overlapSize <- 0
                    totalWidth <-
                        abs(width(localUpNt) - width(localDownNt))
                    jcDist <- 0
                }

                different3UTRoverlap <-
                    jcDist < ntJCsimCutoff &
                    totalWidth - overlapSize > ntCutoff

                # make repport
                localIndex <-
                    which(
                        isoComparison$featureCompared == '3_utr_seq_similarity'
                    )
                isoComparison$isoformsDifferent[localIndex] <-
                    different3UTRoverlap

                # make description
                if (different3UTRoverlap & addDescription) {
                    utr3Gain <-
                        transcriptData$isoform_id[which.max(
                            transcriptData$utr3length
                        )] == upIso

                    if (utr3Gain) {
                        isoComparison$switchConsequence[localIndex] <-
                            '3UTR is longer'
                    } else {
                        isoComparison$switchConsequence[localIndex] <-
                            '3UTR is shorter'
                    }

                }
            }
        }

        if ('NMD_status'                %in% consequencesToAnalyze) {
            if (all(!is.na(transcriptData$PTC))) {
                ptcDifferent <- transcriptData$PTC[1] != transcriptData$PTC[2]

                # make repport
                localIndex <-
                    which(isoComparison$featureCompared == 'NMD_status')
                isoComparison$isoformsDifferent[localIndex] <-
                    ptcDifferent

                # make description
                if (ptcDifferent & addDescription) {
                    upIsSensitive <-
                        transcriptData$PTC[which(
                            transcriptData$isoform_id == upIso
                        )]

                    if (upIsSensitive) {
                        isoComparison$switchConsequence[localIndex] <-
                            'NMD sensitive'
                    } else {
                        isoComparison$switchConsequence[localIndex] <-
                            'NMD insensitive'
                    }
                }
            }
        }

        if ('coding_potential'          %in% consequencesToAnalyze) {
            if (all(!is.na(transcriptData$codingPotential))) {
                cpDifferent <-
                    transcriptData$codingPotential[1] !=
                    transcriptData$codingPotential[2]

                # make repport
                localIndex <-
                    which(isoComparison$featureCompared == 'coding_potential')
                isoComparison$isoformsDifferent[localIndex] <-
                    cpDifferent

                # make description
                if (cpDifferent & addDescription) {
                    upIsCoding <-
                        transcriptData$codingPotential[which(
                            transcriptData$isoform_id == upIso
                        )]

                    if (upIsCoding) {
                        isoComparison$switchConsequence[localIndex] <-
                            'Transcript is coding'
                    } else {
                        isoComparison$switchConsequence[localIndex] <-
                            'Transcript is Noncoding'
                    }
                }
            }
        }

        if ('domains_identified'        %in% consequencesToAnalyze) {
            if (sum(!is.na(transcriptData$orfTransciptLength)) > 0) {
                domianNames <- lapply(domanDataSplit, function(x)
                    x$hmm_name)

                #t1 <- table(domianNames[[1]])
                #t2 <- table(domianNames[[2]])

                t1 <- rle(domianNames[[1]])
                t2 <- rle(domianNames[[2]])

                differentDomains <- !identical(t1, t2)
                # make repport
                localIndex <-
                    which(isoComparison$featureCompared == 'domains_identified')
                isoComparison$isoformsDifferent[localIndex] <-
                    differentDomains

                if (differentDomains & addDescription) {
                    if (sum(t1$lengths) != sum(t2$lengths)) {
                        nrDomains <- sapply(domanDataSplit, nrow)
                        upHasMoreDomains <-
                            names(nrDomains)[which.max(nrDomains)] == upIso

                        if (upHasMoreDomains) {
                            isoComparison$switchConsequence[localIndex] <-
                                'Domain gain'
                        } else {
                            isoComparison$switchConsequence[localIndex] <-
                                'Domain loss'
                        }
                    } else {
                        isoComparison$switchConsequence[localIndex] <-
                            'Domain switch'
                    }
                }
            }
        }

        if ('genomic_domain_position'   %in% consequencesToAnalyze) {
            if (sum(!is.na(transcriptData$orfTransciptLength)) > 0) {
                localDomanDataSplit <-
                    lapply(domanDataSplit, function(aDF) {
                        aDF[
                            sort.list(aDF$pfamStartGenomic),
                            c('hmm_name','pfamStartGenomic','pfamEndGenomic')
                        ]
                    })
                genomicDomainDifferent <-
                    !identical(
                        localDomanDataSplit[[1]],
                        localDomanDataSplit[[2]]
                    )

                # make repport
                localIndex <-
                    which(
                        isoComparison$featureCompared ==
                            'genomic_domain_position'
                    )
                isoComparison$isoformsDifferent[localIndex] <-
                    genomicDomainDifferent
            }
        }

        if ('domain_length'             %in% consequencesToAnalyze) {
            if (sum(!is.na(transcriptData$orfTransciptLength)) > 0) {
                domanDataSplit


                domanDataSplit <- lapply(
                    domanDataSplit,
                    function(x) {
                        x$length <- x$orf_aa_end - x$orf_aa_start + 1
                        return(x)
                    }
                )

                ### Analyze overlap
                nDom <- sapply(domanDataSplit, nrow)
                if( all(nDom) ) {
                    domainRanges <- lapply(
                        domanDataSplit,
                        function(x) {
                            if( x$pfamStartGenomic[1] < x$pfamEndGenomic[1]) {
                                GRanges(
                                    seqnames = 'artificial',
                                    IRanges(
                                        start = x$pfamStartGenomic,
                                        end   = x$pfamEndGenomic
                                    ),
                                    hmm_name  = x$hmm_name,
                                    orfLength = x$length
                                )
                            } else {
                                GRanges(
                                    seqnames = 'artificial',
                                    IRanges(
                                        start = x$pfamEndGenomic,
                                        end   = x$pfamStartGenomic
                                    ),
                                    hmm_name  = x$hmm_name,
                                    orfLength = x$length
                                )
                            }

                        }
                    )

                    ### Calculate overlap
                    localOverlap <- findOverlaps(
                        query   = domainRanges[[ upIso   ]],
                        subject = domainRanges[[ downIso ]]
                    )

                    if(length(localOverlap)) {

                        localOverlapDf <- as.data.frame(localOverlap)
                        localOverlapDf$upName <- domainRanges[[ upIso   ]]$hmm_name[queryHits   (localOverlap)]
                        localOverlapDf$dnName <- domainRanges[[ downIso ]]$hmm_name[subjectHits (localOverlap)]
                        localOverlapDf$upLength <- domainRanges[[ upIso   ]]$orfLength[queryHits   (localOverlap)]
                        localOverlapDf$dnLength <- domainRanges[[ downIso ]]$orfLength[subjectHits (localOverlap)]

                        ### Subset to same domain
                        localOverlapDf <- localOverlapDf[which(
                            localOverlapDf$upName == localOverlapDf$dnName
                        ),]

                        localOverlapDf$lengthDiff <-
                            abs(localOverlapDf$upLength - localOverlapDf$dnLength) > AaCutoff

                        if (!is.null(AaFracCutoff)) {
                            localOverlapDf$minLength <- apply(localOverlapDf[,c('upLength','dnLength')], 1, min)
                            localOverlapDf$maxLength <- apply(localOverlapDf[,c('upLength','dnLength')], 1, max)

                            localOverlapDf$lengthDiff <-
                                localOverlapDf$minLength / localOverlapDf$maxLength < AaFracCutoff & localOverlapDf$lengthDiff

                        }

                        localOverlapDfDiff <- localOverlapDf[which(
                            localOverlapDf$lengthDiff
                        ),]

                        differentDomainLength <- nrow(localOverlapDfDiff) > 0

                        localIndex <-
                            which(isoComparison$featureCompared == 'domain_length')
                        isoComparison$isoformsDifferent[localIndex] <-
                            differentDomainLength

                    } else {
                        differentDomainLength <- FALSE
                    }

                } else {
                    differentDomainLength <- FALSE
                }

                ### Repport if any difference
                if (differentDomainLength & addDescription) {

                    localOverlapDfDiff$maxIsUp <- localOverlapDfDiff$maxLength == localOverlapDfDiff$upLength

                    ### Deside consequence
                    if(
                        all( c(TRUE, FALSE) %in% localOverlapDfDiff$maxIsUp )
                    ) {
                        isoComparison$switchConsequence[localIndex] <-
                            #'IDR length gain and loss'
                            'Mixed Domain length differences'
                    } else if(
                        all(localOverlapDfDiff$maxIsUp)
                    ) {
                        isoComparison$switchConsequence[localIndex] <-
                            'Domain length gain'
                    } else if(
                        all( ! localOverlapDfDiff$maxIsUp )
                    ) {
                        isoComparison$switchConsequence[localIndex] <-
                            'Domain length loss'
                    } else {
                        stop('Something with Domain length analysis went wrong')
                    }

                }

            }
        }

        if ('IDR_identified'            %in% consequencesToAnalyze) {
            if (sum(!is.na(transcriptData$orfTransciptLength)) > 0) {

                nIdr <- sapply(idrDataSplit, nrow)
                if( all(nIdr) ) {
                    idrRanges <- lapply(
                        idrDataSplit,
                        function(x) {
                            if( x$idrStartGenomic[1] < x$idrEndGenomic[1]) {
                                IRanges(
                                    start = x$idrStartGenomic,
                                    end   = x$idrEndGenomic
                                )
                            } else {
                                IRanges(
                                    start = x$idrEndGenomic,
                                    end   = x$idrStartGenomic
                                )
                            }

                        }
                    )

                    ### Calculate overlap
                    overlap1 <- overlapsAny(idrRanges[[1]], idrRanges[[2]])
                    overlap2 <- overlapsAny(idrRanges[[2]], idrRanges[[1]])
                    #overlap1 <- grangesFracOverlap(idrRanges[[1]], idrRanges[[2]])$fracOverlap >= maxIdrFracOverlap
                    #overlap2 <- grangesFracOverlap(idrRanges[[2]], idrRanges[[1]])$fracOverlap >= maxIdrFracOverlap

                } else if( nIdr[1] == 0 & nIdr[2] == 0 ) {
                    overlap1 <- NA
                    overlap2 <- NA
                } else if( nIdr[1] == 0 ) {
                    overlap1 <- TRUE
                    overlap2 <- FALSE
                } else {
                    overlap1 <- FALSE
                    overlap2 <- TRUE
                }

                ### Test overlap
                differentIdr <- any( c(
                    ! overlap1,
                    ! overlap2
                ), na.rm = TRUE)

                # make repport
                localIndex <-
                    which(isoComparison$featureCompared == 'IDR_identified')
                isoComparison$isoformsDifferent[localIndex] <-
                    differentIdr

                if (differentIdr & addDescription) {

                    if( any(!overlap1) & any(!overlap2) ) {
                        isoComparison$switchConsequence[localIndex] <-
                            'IDR switch'
                    } else {
                        if( any(!overlap1) ) {
                            theDiff <- 1
                        } else {
                            theDiff <- 2
                        }

                        upHasMoreIDR <- names(idrDataSplit)[theDiff] == upIso

                        if (upHasMoreIDR) {
                            isoComparison$switchConsequence[localIndex] <-
                                'IDR gain'
                        } else {
                            isoComparison$switchConsequence[localIndex] <-
                                'IDR loss'
                        }
                    }
                }
            }
        }

        if ('IDR_type'                  %in% consequencesToAnalyze) {
            if (sum(!is.na(transcriptData$orfTransciptLength)) > 0) {

                nIdr <- sapply(idrDataSplit, nrow)
                if( all(nIdr) ) {
                    idrRanges <- lapply(
                        idrDataSplit,
                        function(x) {
                            if( x$idrStartGenomic[1] < x$idrEndGenomic[1]) {
                                GRanges(
                                    seqnames = 'artificial',
                                    IRanges(
                                        start = x$idrStartGenomic,
                                        end   = x$idrEndGenomic
                                    ),
                                    type  = x$idr_type
                                )
                            } else {
                                GRanges(
                                    seqnames = 'artificial',
                                    IRanges(
                                        start = x$idrEndGenomic,
                                        end   = x$idrStartGenomic
                                    ),
                                    type  = x$idr_type
                                )
                            }

                        }
                    )

                    ### Calculate overlap
                    localOverlap <- findOverlaps(
                        query   = idrRanges[[ upIso   ]],
                        subject = idrRanges[[ downIso ]]
                    )

                    if(length(localOverlap)) {
                        localOverlapDf <- as.data.frame(localOverlap)
                        localOverlapDf$upType <- idrRanges[[ upIso   ]]$type[queryHits   (localOverlap)]
                        localOverlapDf$dnType <- idrRanges[[ downIso ]]$type[subjectHits (localOverlap)]

                        differentIdrType <- any( na.omit(
                            localOverlapDf$upType != localOverlapDf$dnType
                        ))

                        localIndex <-
                            which(isoComparison$featureCompared == 'IDR_type')
                        isoComparison$isoformsDifferent[localIndex] <-
                            differentIdrType
                    } else {
                        differentIdrType <- FALSE
                    }
                } else {
                    differentIdrType <- FALSE
                }

                if (differentIdrType & addDescription) {

                    localOverlapDf <- localOverlapDf[which(
                        localOverlapDf$upType != localOverlapDf$dnType
                    ),]

                    localOverlapDf$bindingGain <-
                        localOverlapDf$dnType == 'IDR' & localOverlapDf$upType == 'IDR_w_binding_region'

                    localOverlapDf$bindingLoss <-
                        localOverlapDf$dnType == 'IDR_w_binding_region' & localOverlapDf$upType == 'IDR'

                    ### Deside consequence
                    if(
                        any( localOverlapDf$bindingGain) &
                        any(localOverlapDf$bindingLoss)
                    ) {
                        isoComparison$switchConsequence[localIndex] <-
                            'IDR w binding region switch'
                    } else if(
                          any( localOverlapDf$bindingGain ) &
                        ! any( localOverlapDf$bindingLoss )
                    ) {
                        isoComparison$switchConsequence[localIndex] <-
                            'IDR w binding region gain'
                    } else if(
                        ! any( localOverlapDf$bindingGain ) &
                        any( localOverlapDf$bindingLoss )
                    ) {
                        isoComparison$switchConsequence[localIndex] <-
                            'IDR w binding region loss'
                    } else {
                        stop('Something with idr binding regions went wrong')
                    }

                }
            }
        }

        if ('IDR_length'                %in% consequencesToAnalyze) {
            if (sum(!is.na(transcriptData$orfTransciptLength)) > 0) {

                idrDataSplit <- lapply(
                    idrDataSplit,
                    function(x) {
                        x$length <- x$orf_aa_end - x$orf_aa_start + 1
                        return(x)
                    }
                )

                ### Analyze overlap
                nIdr <- sapply(idrDataSplit, nrow)
                if( all(nIdr) ) {
                    idrRanges <- lapply(
                        idrDataSplit,
                        function(x) {
                            if( x$idrStartGenomic[1] < x$idrEndGenomic[1]) {
                                GRanges(
                                    seqnames = 'artificial',
                                    IRanges(
                                        start = x$idrStartGenomic,
                                        end   = x$idrEndGenomic
                                    ),
                                    type  = x$idr_type,
                                    orfLength = x$length
                                )
                            } else {
                                GRanges(
                                    seqnames = 'artificial',
                                    IRanges(
                                        start = x$idrEndGenomic,
                                        end   = x$idrStartGenomic
                                    ),
                                    type  = x$idr_type,
                                    orfLength = x$length
                                )
                            }

                        }
                    )

                    ### Calculate overlap
                    localOverlap <- findOverlaps(
                        query   = idrRanges[[ upIso   ]],
                        subject = idrRanges[[ downIso ]]
                    )

                    if(length(localOverlap)) {

                        localOverlapDf <- as.data.frame(localOverlap)
                        localOverlapDf$upLength <- idrRanges[[ upIso   ]]$orfLength[queryHits   (localOverlap)]
                        localOverlapDf$dnLength <- idrRanges[[ downIso ]]$orfLength[subjectHits (localOverlap)]


                        localOverlapDf$lengthDiff <-
                            abs(localOverlapDf$upLength - localOverlapDf$dnLength) > AaCutoff

                        if (!is.null(AaFracCutoff)) {
                            localOverlapDf$minLength <- apply(localOverlapDf[,c('upLength','dnLength')], 1, min)
                            localOverlapDf$maxLength <- apply(localOverlapDf[,c('upLength','dnLength')], 1, max)

                            localOverlapDf$lengthDiff <-
                                localOverlapDf$minLength / localOverlapDf$maxLength < AaFracCutoff & localOverlapDf$lengthDiff

                        }


                        localOverlapDfDiff <- localOverlapDf[which(
                            localOverlapDf$lengthDiff
                        ),]

                        differentIdrLength <- nrow(localOverlapDfDiff) > 0

                        localIndex <-
                            which(isoComparison$featureCompared == 'IDR_length')
                        isoComparison$isoformsDifferent[localIndex] <-
                            differentIdrLength

                    } else {
                        differentIdrType <- FALSE
                    }

                } else {
                    differentIdrType <- FALSE
                }

                ### Repport if any difference
                if (differentIdrType & addDescription) {

                    localOverlapDfDiff$maxIsUp <- localOverlapDfDiff$maxLength == localOverlapDfDiff$upLength

                    ### Deside consequence
                    if(
                        all( c(TRUE, FALSE) %in% localOverlapDfDiff$maxIsUp )
                    ) {
                        isoComparison$switchConsequence[localIndex] <-
                            #'IDR length gain and loss'
                            'Mixed IDR length differences'
                    } else if(
                        all(localOverlapDfDiff$maxIsUp)
                    ) {
                        isoComparison$switchConsequence[localIndex] <-
                            'IDR length gain'
                    } else if(
                        all( ! localOverlapDfDiff$maxIsUp )
                    ) {
                        isoComparison$switchConsequence[localIndex] <-
                            'IDR length loss'
                    } else {
                        stop('Something with idr length went wrong')
                    }

                }
            }
        }

        if ('signal_peptide_identified' %in% consequencesToAnalyze) {
            if (sum(!is.na(transcriptData$orfTransciptLength)) > 0) {
                nrSignalPeptides <- sapply(peptideDataSplit, nrow)

                differentSignaPeptied <-
                    nrSignalPeptides[1] != nrSignalPeptides[2]

                # make repport
                localIndex <-
                    which(
                        isoComparison$featureCompared ==
                            'signal_peptide_identified'
                    )
                isoComparison$isoformsDifferent[localIndex] <-
                    differentSignaPeptied

                if (differentSignaPeptied & addDescription) {
                    upHasSignal <-
                        names(nrSignalPeptides)[which.max(
                            nrSignalPeptides
                        )] == upIso

                    if (upHasSignal) {
                        isoComparison$switchConsequence[localIndex] <-
                            'Signal peptide gain'
                    } else {
                        isoComparison$switchConsequence[localIndex] <-
                            'Signal peptide loss'
                    }
                }
            }
        }

        if ('sub_cell_location'         %in% consequencesToAnalyze) {
            if (sum(!is.na(transcriptData$orfTransciptLength)) == 2) {
                if( sum(!is.na(transcriptData$sub_cell_location)) == 2 ) {
                    differentLoc <-
                        transcriptData$sub_cell_location[1] != transcriptData$sub_cell_location[2]

                    # make repport
                    localIndex <-
                        which(
                            isoComparison$featureCompared ==
                                'sub_cell_location'
                        )
                    isoComparison$isoformsDifferent[localIndex] <-
                        differentLoc

                    if (differentLoc & addDescription) {
                        upLoc <- transcriptData$sub_cell_location[which(
                            transcriptData$isoform_id == upIso
                        )]
                        dnLoc <- transcriptData$sub_cell_location[which(
                            transcriptData$isoform_id != upIso
                        )]

                        #isoComparison$switchConsequence[localIndex] <- paste0('Location switch to ', upLoc)
                        isoComparison$switchConsequence[localIndex] <- paste0(
                            'Location switch from ',
                            dnLoc,
                            ' to ', upLoc
                        )
                    }
                }
            }
        }

        if ('solubility_status'         %in% consequencesToAnalyze) {
            if (sum(!is.na(transcriptData$orfTransciptLength)) == 2) {
                if( sum(!is.na(transcriptData$solubility_status)) == 2 ) {
                    differentSolubility <-
                        transcriptData$solubility_status[1] !=
                        transcriptData$solubility_status[2]

                    # make repport
                    localIndex <-
                        which(
                            isoComparison$featureCompared ==
                                'solubility_status'
                        )
                    isoComparison$isoformsDifferent[localIndex] <-
                        differentSolubility

                    if (differentSolubility & addDescription) {
                        upSoluble <-
                            'Soluble' ==
                            transcriptData$solubility_status[which(
                                transcriptData$isoform_id == upIso
                            )]
                        if (upSoluble) {
                            isoComparison$switchConsequence[localIndex] <-
                                'Membrane tethering loss'
                        } else {
                            isoComparison$switchConsequence[localIndex] <-
                                'Membrane tethering gain'
                        }
                    }
                }
            }

        }

    }

    if (onlyRepportDifferent) {
        isoComparison <-
            isoComparison[which(
                isoComparison$isoformsDifferent
        ),]
    }

    return(isoComparison)
}


### Summarizing consequences
extractConsequenceSummary <- function(
    switchAnalyzeRlist,
    consequencesToAnalyze = 'all',
    includeCombined = FALSE,
    asFractionTotal = FALSE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    plot = TRUE,
    plotGenes = FALSE,
    simplifyLocation = TRUE,
    localTheme = theme_bw(),
    returnResult = FALSE
) {
    ### Test input
    if (TRUE) {
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(
                'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
            )
        }

        # test wether switching have been analyzed
        if (!any(!is.na(
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value
        ))) {
            stop(
                'The analsis of isoform switching must be performed before functional consequences can be analyzed. Please run detectIsoformSwitching() and try again.'
            )
        }
        # test whether switches have been predicted
        if (is.null(switchAnalyzeRlist$switchConsequence)) {
            stop(
                'The analsis of isoform switch consequences must be performed before it can be summarized. Please use analyzeSwitchConsequences() and try again.'
            )
        }

        # input format
        if (dIFcutoff < 0 | dIFcutoff > 1) {
            stop('The dIFcutoff must be in the interval [0,1].')
        }
        if (alpha < 0 |
            alpha > 1) {
            warning('The alpha parameter should usually be between 0 and 1 ([0,1]).')
        }
        if (alpha > 0.05) {
            warning(
                'Most journals and scientists consider an alpha larger than 0.05 untrustworthy. We therefore recommend using alpha values smaller than or queal to 0.05'
            )
        }

        ### Consequences to analyze
        acceptedTypes <- c(
            # Transcript
            'tss',
            'tts',
            'last_exon',
            'isoform_length',
            'exon_number',
            'intron_structure',
            'intron_retention',
            'isoform_class_code',
            # cpat
            'coding_potential',
            # ORF
            'ORF_genomic',
            'ORF_length',
            '5_utr_length',
            '3_utr_length',
            # seq similarity
            'isoform_seq_similarity',
            'ORF_seq_similarity',
            '5_utr_seq_similarity',
            '3_utr_seq_similarity',
            # ORF
            'NMD_status',
            # pfam
            'domains_identified',
            'genomic_domain_position',
            'domain_length',

            # SignalIP
            'signal_peptide_identified',

            # IDR
            'IDR_identified',
            'IDR_length',
            'IDR_type',

            # sub cell
            'sub_cell_location',
            'solubility_status'
        )

        if (!all(consequencesToAnalyze %in% c('all', acceptedTypes))) {
            stop(
                'The argument(s) supplied to \'typeOfconsequence\' are not accepted. Please see ?analyzeSwitchConsequences under details for description of which strings are allowed.'
            )
        }

        consequencesAnalyzed <-
            unique(switchAnalyzeRlist$switchConsequence$featureCompared)
        if ('all' %in% consequencesToAnalyze) {
            consequencesToAnalyze <- consequencesAnalyzed
        }

        consequencesNotAnalyzed <-
            setdiff(consequencesToAnalyze, consequencesAnalyzed)
        if (length(consequencesNotAnalyzed)) {
            warning(
                paste(
                    'The following consequences appear not to have been analyzed and will therefor not be summarized:',
                    paste(consequencesNotAnalyzed, collapse = ', '),
                    sep = ' '
                )
            )
        }


    }

    ### Massage for plotting
    localSwitchConsequences <-
        switchAnalyzeRlist$switchConsequence[which(
            switchAnalyzeRlist$switchConsequence$isoformsDifferent &
                switchAnalyzeRlist$switchConsequence$featureCompared %in%
                consequencesToAnalyze
        ),]
    if (!nrow(localSwitchConsequences)) {
        stop('No swithces with consequences were found')
    }

    localSwitchConsequences <-
        localSwitchConsequences[which(!is.na(
            localSwitchConsequences$switchConsequence
        )), ]
    if (!nrow(localSwitchConsequences)) {
        stop('No swithces with describable consequences were found')
    }
    localSwitchConsequences$Comparison <-
        paste(
            localSwitchConsequences$condition_1,
            'vs',
            localSwitchConsequences$condition_2,
            sep = ' '
        )

    if (includeCombined) {
        tmp <- localSwitchConsequences

        tmp$featureCompared <- 'combined'
        tmp$switchConsequence <- 'any consequence'
        tmp <- unique(tmp)

        localSwitchConsequences <-
            rbind(localSwitchConsequences, tmp)
    }

    ### Handle location
    if( simplifyLocation & 'sub_cell_location' %in% localSwitchConsequences$featureCompared ) {
        toModifyIndex <- which(localSwitchConsequences$featureCompared == 'sub_cell_location')

        localSwitchConsequences$switchConsequence[toModifyIndex] <- sapply(
            strsplit(
                localSwitchConsequences$switchConsequence[toModifyIndex],
                ' from | to '
            ),
            function(x) {
                paste(x[c(1,3)], collapse = ' to ')
            }
        )
    }


    ### Extract Sig iso
    isoResTest <-
        any(!is.na(
            switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
        ))
    if (isoResTest) {
        sigIso <- switchAnalyzeRlist$isoformFeatures[which(
            switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value < alpha &
                abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
        ),
        c('iso_ref', 'gene_ref')]
    } else {
        sigIso <- switchAnalyzeRlist$isoformFeatures[which(
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value < alpha &
                abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
        ),
        c('iso_ref', 'gene_ref')]
    }

    ### summarize count
    myNumbers <-
        plyr::ddply(
            localSwitchConsequences,
            .variables = c(
                'switchConsequence',
                'featureCompared',
                'Comparison'
            ),
            .drop = TRUE,
            .fun = function(aDF) {
                data.frame(
                    featureCompared = aDF$featureCompared[1],
                    nrGenesWithConsequences = length(
                        intersect(aDF$gene_ref, sigIso$gene_ref)
                    ),
                    nrIsoWithConsequences = length(intersect(
                        c(aDF$iso_ref_up, aDF$iso_ref_down), sigIso$iso_ref
                    )),
                    stringsAsFactors = FALSE
                )
            }
        )

    myNumbers$featureCompared <-
        startCapitalLetter(gsub('_', ' ', myNumbers$featureCompared))
    myNumbers <-
        myNumbers[, c(
            'Comparison',
            'featureCompared',
            'switchConsequence',
            'nrGenesWithConsequences',
            'nrIsoWithConsequences'
        )]

    if (asFractionTotal) {
        totalNrSwitches <-
            extractSwitchSummary(switchAnalyzeRlist,
                                 alpha = alpha,
                                 dIFcutoff = dIFcutoff)

        ### add to numbers
        matchVec <-
            match(myNumbers$Comparison, totalNrSwitches$Comparison)
        myNumbers$nrIsoforms <- totalNrSwitches$nrIsoforms[matchVec]
        myNumbers$nrGenes <- totalNrSwitches$nrGenes[matchVec]

        myNumbers$geneFraction <-
            round(myNumbers$nrGenesWithConsequences / myNumbers$nrGenes,
                  digits = 4)
        myNumbers$isoFraction  <-
            round(myNumbers$nrIsoWithConsequences   / myNumbers$nrIsoforms,
                  digits = 4)
    }

    if (plot) {
        myNumbers$plotComparison <-
            gsub(' vs ', '\nvs\n', myNumbers$Comparison)

        myNumbers$featureCompared <- newLineAtMiddel(
            myNumbers$featureCompared
        )

        ### Make plot
        if (plotGenes) {
            if (asFractionTotal) {
                g1 <- ggplot(myNumbers, aes(
                    x = switchConsequence,
                    y = geneFraction
                )) +
                    labs(
                        x = 'Consequence of switch\n(features of the upregulated isoform)',
                        y = 'Fraction of genes'
                    )
            } else {
                g1 <-
                    ggplot(myNumbers, aes(
                        x = switchConsequence,
                        y = nrGenesWithConsequences
                    )) +
                    labs(
                        x = 'Consequence of switch\n(features of the upregulated isoform)',
                        y = 'Number of genes'
                        )
            }
        } else {
            if (asFractionTotal) {
                g1 <- ggplot(
                    myNumbers,
                    aes(x = switchConsequence, y = isoFraction)
                ) + labs(
                    x = 'Consequence of switch\n(features of the upregulated isoform)',
                    y = 'Fraction of isoforms')
            } else {
                g1 <- ggplot(
                    myNumbers,
                    aes(x = switchConsequence, y = nrIsoWithConsequences)
                ) + labs(
                    x = 'Consequence of switch\n(features of the upregulated isoform)',
                    y = 'Number of isoforms'
                )
            }
        }

        g1 <- g1 + geom_bar(stat = "identity") +
            facet_grid(plotComparison ~ featureCompared,
                       scales = 'free',
                       space = 'free_x') +
            localTheme +
            theme(strip.text.y = element_text(angle = 0)) +
            theme(axis.text.x = element_text(
                angle = -45,
                hjust = 0,
                vjust = 1
            ))

        print(g1)

        myNumbers$plotComparison <- NULL
    }

    if (returnResult) {
        return(myNumbers)
    }

}

extractConsequenceEnrichment <- function(
    switchAnalyzeRlist,
    consequencesToAnalyze = 'all',
    alpha=0.05,
    dIFcutoff = 0.1,
    countGenes = TRUE,
    analysisOppositeConsequence=FALSE,
    plot=TRUE,
    localTheme = theme_bw(base_size = 12),
    minEventsForPlotting = 10,
    returnResult=TRUE,
    returnSummary=TRUE
) {
    ### Test input
    if (TRUE) {
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(
                'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
            )
        }

        # test wether switching have been analyzed
        if (!any(!is.na(
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value
        ))) {
            stop(
                'The analsis of isoform switching must be performed before functional consequences can be analyzed. Please run detectIsoformSwitching() and try again.'
            )
        }
        # test whether switches have been predicted
        if (is.null(switchAnalyzeRlist$switchConsequence)) {
            stop(
                'The analsis of isoform switch consequences must be performed before it can be summarized. Please use analyzeSwitchConsequences() and try again.'
            )
        }

        # input format
        if (dIFcutoff < 0 | dIFcutoff > 1) {
            stop('The dIFcutoff must be in the interval [0,1].')
        }
        if (alpha < 0 |
            alpha > 1) {
            warning('The alpha parameter should usually be between 0 and 1 ([0,1]).')
        }
        if (alpha > 0.05) {
            warning(
                'Most journals and scientists consider an alpha larger than 0.05 untrustworthy. We therefore recommend using alpha values smaller than or queal to 0.05'
            )
        }

        ### Consequences to analyze
        acceptedTypes <- c(
            # Transcript
            'tss',
            'tts',
            'last_exon',
            'isoform_length',
            'exon_number',
            'intron_structure',
            'intron_retention',
            'isoform_class_code',
            # cpat
            'coding_potential',
            # ORF
            'ORF_genomic',
            'ORF_length',
            '5_utr_length',
            '3_utr_length',
            # seq similarity
            'isoform_seq_similarity',
            'ORF_seq_similarity',
            '5_utr_seq_similarity',
            '3_utr_seq_similarity',
            # ORF
            'NMD_status',
            # pfam
            'domains_identified',
            'genomic_domain_position',
            'domain_length',

            # SignalIP
            'signal_peptide_identified',

            # IDR
            'IDR_identified',
            'IDR_length',
            'IDR_type',

            # sub cell
            'sub_cell_location',
            'solubility_status'
        )

        if (!all(consequencesToAnalyze %in% c('all', acceptedTypes))) {
            stop(
                'The argument(s) supplied to \'typeOfconsequence\' are not accepted. Please see ?summarizeSwitchConsequences under details for description of which strings are allowed.'
            )
        }

        consequencesAnalyzed <-
            unique(switchAnalyzeRlist$switchConsequence$featureCompared)
        if ('all' %in% consequencesToAnalyze) {
            consequencesToAnalyze <- consequencesAnalyzed
        }

        consequencesNotAnalyzed <-
            setdiff(consequencesToAnalyze, consequencesAnalyzed)
        if (length(consequencesNotAnalyzed)) {
            warning(
                paste(
                    'The following consequences appear not to have been analyzed and will therefor not be summarized:',
                    paste(consequencesNotAnalyzed, collapse = ', '),
                    sep = ' '
                )
            )
        }

        hasLocation <- 'sub_cell_location' %in%
            colnames(switchAnalyzeRlist$isoformFeatures)


    }

    ### Extract non-location consequences
    if(TRUE) {
        localConseq <- switchAnalyzeRlist$switchConsequence[
            which( !is.na(
                switchAnalyzeRlist$switchConsequence$switchConsequence
            ))
            ,]
        localConseq <- localConseq[which(
            ! grepl('switch', localConseq$switchConsequence)
        ),]

        ### make list with levels
        levelList <- list(
            tss=c('Tss more upstream','Tss more downstream'),
            tts=c('Tts more downstream','Tts more upstream'),
            last_exon=c('Last exon more downstream','Last exon more upstream'),
            isoform_length=c('Length gain','Length loss'),
            isoform_seq_similarity=c('Length gain','Length loss'),
            exon_number=c('Exon gain','Exon loss'),
            intron_retention=c('Intron retention gain','Intron retention loss'),
            ORF_length=c('ORF is longer','ORF is shorter'),
            ORF=c('Complete ORF loss','Complete ORF gain'),
            x5_utr_length=c('5UTR is longer','5UTR is shorter'),
            x3_utr_length=c('3UTR is longer','3UTR is shorter'),
            NMD_status=c('NMD sensitive','NMD insensitive'),
            coding_potential=c('Transcript is coding','Transcript is Noncoding'),
            domains_identified=c('Domain gain','Domain loss'),
            domain_length=c('Domain length gain','Domain length loss'),
            IDR_identified = c('IDR gain','IDR loss'),
            IDR_length = c('IDR length gain', 'IDR length loss'),
            IDR_type = c('IDR w binding region gain', 'IDR w binding region loss'),
            signal_peptide_identified=c('Signal peptide gain','Signal peptide loss'),
            solubility_status = c('Membrane tethering gain','Membrane tethering loss')
        )
        levelListDf <- plyr::ldply(levelList, function(x) data.frame(feature=x, stringsAsFactors = FALSE))

        ### Add consequence paris
        localConseq$conseqPair <- levelListDf$.id[match(localConseq$switchConsequence, levelListDf$feature)]
        localConseq <- localConseq[which( !is.na(localConseq$conseqPair)),]

        ### Subset to consequences analyzed
        localConseq <- localConseq[which(
            localConseq$featureCompared %in% consequencesToAnalyze
        ),]
    }

    ### Extract location consequecnes
    if(hasLocation) {
        localLocConseq <- switchAnalyzeRlist$switchConsequence[
            which( !is.na(
                switchAnalyzeRlist$switchConsequence$switchConsequence
            ))
            ,]
        localLocConseq <- localLocConseq[which(
            grepl('Location switch', localLocConseq$switchConsequence)
        ),]

        locationsList <- lapply(
            strsplit(
                localLocConseq$switchConsequence,
                ' from | to '
            ),
            function(x) {
                x[2:3]
            }
        )
        localLocConseq$from <- sapply(
            locationsList,
            function(x) x[1]
        )
        localLocConseq$to <- sapply(
            locationsList,
            function(x) x[2]
        )

    }


    ### Subset to significant features
    if(TRUE) {
        ### Extract Sig iso
        isoResTest <-
            any(!is.na(
                switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
            ))
        if (isoResTest) {
            sigIso <- switchAnalyzeRlist$isoformFeatures[which(
                switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value < alpha &
                    abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
            ),
            c('iso_ref', 'gene_ref')]
        } else {
            sigIso <- switchAnalyzeRlist$isoformFeatures[which(
                switchAnalyzeRlist$isoformFeatures$gene_switch_q_value < alpha &
                    abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
            ),
            c('iso_ref', 'gene_ref')]
        }

        if(isoResTest) {
            localConseq <- localConseq[which(
                localConseq$iso_ref_down %in% sigIso$iso_ref |
                localConseq$iso_ref_up   %in% sigIso$iso_ref
            ),]

            if(hasLocation) {
                localLocConseq <- localLocConseq[which(
                    localLocConseq$iso_ref_down %in% sigIso$iso_ref |
                        localLocConseq$iso_ref_up   %in% sigIso$iso_ref
                ),]
            }
        } else {
            localConseq <- localConseq[which(
                localConseq$gene_ref %in% sigIso$gene_ref
            ),]

            if(hasLocation) {
                localLocConseq <- localLocConseq[which(
                    localLocConseq$gene_ref %in% sigIso$gene_ref
                ),]
            }

        }
    }


    ### Summarize gain vs loss for each consequence in each condition
    consequenceBalance <- plyr::ddply(
        .data = localConseq,
        .variables = c('condition_1','condition_2','conseqPair'),
        #.inform = TRUE,
        .fun = function(aDF) { # aDF <- localConseq[1:20,]
            ### Add levels
            if(analysisOppositeConsequence) {
                localLvl <- rev(sort(
                    levelList[[ aDF$conseqPair[1] ]]
                ))
            } else {
                localLvl <- sort(
                    levelList[[ aDF$conseqPair[1] ]]
                )
            }
            aDF$switchConsequence <- factor(
                aDF$switchConsequence,
                levels=localLvl
            )

            ### Summarize category
            if( countGenes ) {
                df2 <- aDF[
                    which(!is.na(aDF$switchConsequence)),
                    c('gene_id','switchConsequence')
                ]
                localNumber <- plyr::ddply(df2, .drop = FALSE, .variables = 'switchConsequence', function(x) {
                    data.frame(
                        Freq = length(unique(x$gene_id))
                    )
                })
                colnames(localNumber)[1] <- 'Var1'
            } else {
                localNumber <- as.data.frame(table(aDF$switchConsequence))
            }

            if(nrow(localNumber) == 2) {
                localTest <- suppressWarnings(
                    prop.test(localNumber$Freq[1], sum(localNumber$Freq))
                )

                localRes <- data.frame(
                    feature=stringr::str_c(
                        localNumber$Var1[1],
                        ' (paired with ',
                        localNumber$Var1[2],
                        ')'
                    ),
                    propOfRelevantEvents=localTest$estimate,
                    stringsAsFactors = FALSE
                )

                localRes$propCiLo <- min(localTest$conf.int)
                localRes$propCiHi <- max(localTest$conf.int)
                localRes$propPval <- localTest$p.value
            } else {
                warning('Somthing strange happend - contact developer with reproducible example')
            }

            localRes$nUp   <- localNumber$Freq[which( localNumber$Var1 == levels(localNumber$Var1)[1] )]
            localRes$nDown <- localNumber$Freq[which( localNumber$Var1 == levels(localNumber$Var1)[2] )]

            return(localRes)
        }
    )

    ### Summarize gain vs loss for each location
    if( hasLocation ) {
        locationBalance <- plyr::ddply(
            .data = localLocConseq,
            .variables = c('condition_1','condition_2'),
            #.inform = TRUE,
            .fun = function(aDF) { # aDF <- localLocConseq[1:20,]
                locationsAnalyzed <- unique(c(
                    aDF$from, aDF$to
                ))
                locationsAnalyzedList <- split(
                    locationsAnalyzed,
                    locationsAnalyzed
                )

                localCount <- plyr::ldply(
                    locationsAnalyzedList,
                    function(aLocation) {
                        ### Summarize category
                        if( countGenes ) {
                            localNumber <- data.frame(
                                Var1 = paste(c('Location switch to','Location switch away from'), aLocation),
                                Freq = c(
                                    length(unique( aDF$gene_ref[which(aDF$to   == aLocation)])),
                                    length(unique( aDF$gene_ref[which(aDF$from == aLocation)]))
                                )
                            )
                        } else {
                            localNumber <- data.frame(
                                Var1 = paste(c('Location switch to','Location switch away from'), aLocation),
                                Freq = c(
                                    sum(aDF$to   == aLocation),
                                    sum(aDF$from == aLocation)
                                )
                            )
                        }

                        if(nrow(localNumber) == 2) {
                            localTest <- suppressWarnings(
                                prop.test(localNumber$Freq[1], sum(localNumber$Freq))
                            )

                            localRes <- data.frame(
                                feature=stringr::str_c(
                                    localNumber$Var1[1],
                                    ' (paired with ',
                                    localNumber$Var1[2],
                                    ')'
                                ),
                                propOfRelevantEvents=localTest$estimate,
                                stringsAsFactors = FALSE
                            )

                            localRes$propCiLo <- min(localTest$conf.int)
                            localRes$propCiHi <- max(localTest$conf.int)
                            localRes$propPval <- localTest$p.value
                        } else {
                            warning('Somthing strange happend - contact developer with reproducible example')
                        }

                        localRes$nUp   <- localNumber$Freq[1] # order is always fixed
                        localRes$nDown <- localNumber$Freq[2] # order is always fixed

                        return(localRes)

                    }
                )

                return(localCount)
            }
        )

        colnames(locationBalance)[3] <- c('conseqPair')


        consequenceBalance <- rbind(
            consequenceBalance,
            locationBalance
        )
    }

    consequenceBalance$propQval <- p.adjust(consequenceBalance$propPval, method = 'fdr')
    consequenceBalance$Significant <- consequenceBalance$propQval < alpha
    consequenceBalance$Significant <- factor(
        consequenceBalance$Significant,
        levels=c(FALSE,TRUE)
    )

    ### Plot result
    if(plot) {
        if(countGenes) {
            xText <- 'Fraction of Genes Having the Consequence Indicated\n(of Switches Affected by Either of Opposing Consequences)\n(With 95% Confidence Interval)'
        } else {
            xText <- 'Fraction of Switches Having the Consequence Indicated\n(of Switches Affected by Either of Opposing Consequences)\n(With 95% Confidence Interval)'
        }

        consequenceBalance2 <- consequenceBalance[which(
            (consequenceBalance$nUp + consequenceBalance$nDown) >= minEventsForPlotting
        ),]
        if(nrow(consequenceBalance2) == 0) {
            stop('No features left for ploting after filtering with via "minEventsForPlotting" argument.')
        }

        consequenceBalance2$nTot <- consequenceBalance2$nDown + consequenceBalance2$nUp

        ### Add comparison
        consequenceBalance2$Comparison <- paste(
            consequenceBalance2$condition_1,
            'vs',
            consequenceBalance2$condition_2,
            sep='\n'
        )

        ### Massage
        consequenceBalance2$feature2 <- gsub(' \\(', '\n(', consequenceBalance2$feature)

        consequenceBalance2$feature2 <- factor(
            consequenceBalance2$feature2,
            levels = rev(sort(unique(as.character(consequenceBalance2$feature2))))
        )

        g1 <- ggplot(data=consequenceBalance2, aes(y=feature2, x=propOfRelevantEvents, color=Significant)) +
            #geom_point(size=4) +
            geom_errorbarh(aes(xmax = propCiLo, xmin=propCiHi), height = .3) +
            geom_point(aes(size=nTot)) +
            facet_wrap(~Comparison) +
            geom_vline(xintercept=0.5, linetype='dashed') +
            labs(
                x=xText,
                y='Consequence of Isoform Switch\n(and the opposing consequence)'
            ) +
            localTheme +
            theme(axis.text.x=element_text(angle=-45, hjust = 0, vjust=1)) +
            scale_color_manual(name = paste0('FDR < ', alpha), values=c('black','red'), drop=FALSE) +
            guides(
                color = guide_legend(order=1),
                size = guide_legend(order=2)
            ) +
            coord_cartesian(xlim=c(0,1))

        if( countGenes ) {
            g1 <- g1 + scale_size_continuous(name = 'Genes')
        } else {
            g1 <- g1 + scale_size_continuous(name = 'Switches')
        }

        print(g1)
    }

    if(returnResult) {
        if(returnSummary) {
            consequenceBalance$Comparison <- NULL

            return(consequenceBalance)

        } else {
            return(localConseq)
        }
    }
}

extractConsequenceEnrichmentComparison <- function(
    switchAnalyzeRlist,
    consequencesToAnalyze = 'all',
    alpha=0.05,
    dIFcutoff = 0.1,
    countGenes = TRUE,
    analysisOppositeConsequence=FALSE,
    plot=TRUE,
    localTheme = theme_bw(base_size = 14),
    minEventsForPlotting = 10,
    returnResult=TRUE
) {
    ### Test
    if(TRUE) {
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(
                'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
            )
        }
        nComp <- nrow( unique(
            switchAnalyzeRlist$isoformFeatures[,c('condition_1','condition_2')]
        ))

        if( nComp == 1) {
            stop('Cannot do a contrast of different comparisons since only one comparison is analyzed in the switchAnalyzeRlist.')
        }
    }


    ### Extract splicing enrichment
    conseqCount <- extractConsequenceEnrichment(
        switchAnalyzeRlist = switchAnalyzeRlist,
        consequencesToAnalyze = consequencesToAnalyze,
        alpha = alpha,
        dIFcutoff = dIFcutoff,
        countGenes = countGenes,
        minEventsForPlotting = 0,
        analysisOppositeConsequence=analysisOppositeConsequence,
        plot = FALSE,
        returnResult = TRUE
    )
    conseqCount$Comparison <- stringr::str_c(
        conseqCount$condition_1,
        '\nvs\n',
        conseqCount$condition_2
    )
    conseqCount$nTot <- conseqCount$nDown + conseqCount$nUp

    ### Extract pairs
    conseqPairs <- split(as.character(unique(conseqCount$feature)), unique(conseqCount$feature))

    ### Make each pairwise comparison
    myComparisons <- allPairwiseFeatures( unique(conseqCount$Comparison) )

    ### Loop over each pariwise comparison
    fisherRes <- plyr::ddply(myComparisons, c('var1','var2'), function(localComparison) { # localComparison <- myComparisons[1,]
        ### Loop over each
        localAsRes <- plyr::ldply(conseqPairs, .inform = TRUE, function(localConseq) { # localConseq <- 'Membrane tethering gain (paired with Membrane tethering loss)'
            ### Extract local data
            localSpliceCount1 <- conseqCount[which(
                conseqCount$Comparison %in% c(localComparison$var1, localComparison$var2) &
                    conseqCount$feature == localConseq
            ),]

            if(
                nrow(localSpliceCount1) != 2 |
                any(is.na(localSpliceCount1[,c('nUp','nDown')]))
            ) {
                return(NULL)
            }

            ### Test difference
            fisherResult <- fisher.test(localSpliceCount1[,c('nUp','nDown')])

            ###
            fisherTestResult <- data.frame(odds_ratio=fisherResult$estimate, p_value=fisherResult$p.value, lowCI=fisherResult$conf.int[1], highCI=fisherResult$conf.int[2])
            rownames(fisherTestResult) <- NULL

            localSpliceCount1$fisherOddsRatio <- fisherTestResult$odds_ratio
            localSpliceCount1$fisherPvalue    <- fisherTestResult$p_value

            localSpliceCount1$pair <- 1:2
            localSpliceCount1$forPlotting <- any(
                (localSpliceCount1$nUp + localSpliceCount1$nDown) >= minEventsForPlotting
            )

            return(
                localSpliceCount1[,c('Comparison','propOfRelevantEvents','propCiLo','propCiHi','fisherOddsRatio','fisherPvalue','pair','forPlotting','nTot')]
            )
        })

        colnames(localAsRes)[1] <- 'consequence'

        return(localAsRes)
    })

    ### Add comparison
    fisherRes$comp <- paste(
        gsub('\\n',' ',fisherRes$var1),
        gsub('\\n',' ',fisherRes$var2),
        sep='\ncompared to\n'
    )

    ### Multiple test correction
    # perform correction for pair = 1 (to avoid correcting twice for the same pair)
    tmp <- fisherRes[which(fisherRes$pair == 1),]
    tmp$fisherQvalue <- p.adjust(tmp$fisherPvalue, method = 'fdr')

    fisherRes$fisherQvalue <- tmp$fisherQvalue[match(
        stringr::str_c(fisherRes$var1, fisherRes$var2, fisherRes$consequence),
        stringr::str_c(tmp$var1, tmp$var2, tmp$consequence)
    )]
    fisherRes$Significant <- fisherRes$fisherQvalue < alpha
    fisherRes$Significant <- factor(
        fisherRes$Significant,
        levels=c(FALSE,TRUE)
    )

    ### Plot
    if(plot) {
        fisherRes2 <- fisherRes[which(fisherRes$forPlotting),]
        if(nrow(fisherRes2) == 0) {
            stop('No features left to plot after subsetting with \'minEventsForPlotting\'.')
        }

        fisherRes2$consequence <- gsub(' \\(paired','\n\\(paired', fisherRes2$consequence)
        fisherRes2$consequence <- gsub('with ','with\n', fisherRes2$consequence)

        if(countGenes) {
            xText <- 'Fraction of genes having the consequence indicated\n(of the switches affected by either of opposing consequences)\n(with 95% confidence interval)'
        } else {
            xText <- 'Fraction of switches having the consequence indicated\n(of the switches affected by either of opposing consequences)\n(with 95% confidence interval)'
        }

        g1 <- ggplot(data=fisherRes2, aes(y=Comparison, x=propOfRelevantEvents, color=Significant)) +
            #geom_point(aes(size=4)) +
            geom_point(aes(size=nTot)) +
            geom_errorbarh(aes(xmax = propCiHi, xmin=propCiLo), height = .3) +
            facet_grid(comp~consequence, scales = 'free_y') +
            geom_vline(xintercept=0.5, linetype='dashed') +
            labs(x=xText, y='Comparison') +
            #scale_color_manual('Fraction in\nComparisons\nSignifcantly different', values=c('black','red'), drop=FALSE) +
            scale_color_manual(name = paste0('FDR across\ncomparisons\n< ', alpha), values=c('black','red'), drop=FALSE) +
            guides(
                color = guide_legend(order=1),
                size = guide_legend(order=2)
            ) +
            localTheme +
            theme(axis.text.x=element_text(angle=-45, hjust = 0, vjust=1), strip.text.y = element_text(angle = 0)) +
            coord_cartesian(xlim=c(0,1))

        if( countGenes ) {
            g1 <- g1 + scale_size_continuous(name = 'Genes')
        } else {
            g1 <- g1 + scale_size_continuous(name = 'Switches')
        }

        print(g1)

    }

    if(returnResult) {
        fisherRes$nTot <- NULL

        fisherRes$pair <- stringr::str_c('propUp_comparison_', fisherRes$pair)
        colnames(fisherRes)[which(colnames(fisherRes) == 'propOfRelevantEvents')] <- 'propUp'

        fisherRes2 <- reshape2::dcast(data = fisherRes, comp + consequence ~ pair, value.var=c('propUp'))

        matchIndex <- match(
            stringr::str_c(fisherRes2$comp, fisherRes2$consequence),
            stringr::str_c(tmp$comp, tmp$consequence)
        )
        #fisherRes2$fisherOddsRatio <- tmp$fisherOddsRatio[matchIndex]
        fisherRes2$fisherQvalue    <- tmp$fisherQvalue   [matchIndex]

        fisherRes2$Significant <- fisherRes2$fisherQvalue < alpha
        colnames(fisherRes2)[1] <- 'comparisonsCompared'
        fisherRes2$comparisonsCompared <- gsub('\\n', ' ', fisherRes2$comparisonsCompared)
        return(fisherRes2)
    }

}

extractConsequenceGenomeWide <- function(
    switchAnalyzeRlist,
    featureToExtract = 'isoformUsage',
    annotationToAnalyze = 'all',
    alpha=0.05,
    dIFcutoff = 0.1,
    log2FCcutoff = 1,
    violinPlot=TRUE,
    alphas=c(0.05, 0.001),
    localTheme=theme_bw(),
    plot=TRUE,
    returnResult=TRUE
) {
    ### Test input
    if (TRUE) {
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(
                'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
            )
        }
        if( switchAnalyzeRlist$sourceId == 'preDefinedSwitches' ) {
            stop(
                paste(
                    'The switchAnalyzeRlist is made from pre-defined isoform switches',
                    'which means it is made without defining conditions (as it should be).',
                    '\nThis also means it cannot used to plot conditional expression -',
                    'if that is your intention you need to create a new',
                    'switchAnalyzeRlist with the importRdata() function and start over.',
                    sep = ' '
                )
            )
        }

        if (alpha < 0 |
            alpha > 1) {
            warning('The alpha parameter should usually be between 0 and 1 ([0,1]).')
        }
        if (alpha > 0.05) {
            warning(
                'Most journals and scientists consider an alpha larger than 0.05 untrustworthy. We therefore recommend using alpha values smaller than or queal to 0.05'
            )
        }
        if (dIFcutoff < 0 | dIFcutoff > 1) {
            stop('The dIFcutoff must be in the interval [0,1].')
        }
        if (log2FCcutoff < 0) {
            stop(
                'The log2FCcutoff cannot be negative (as the cutoff is applied to absolute values)'
            )
        }
        if (length(alphas) != 2) {
            stop('A vector of length 2 must be given to the argument \'alphas\'')
        }
        if (any(alphas < 0) |
            any(alphas > 1)) {
            stop('The \'alphas\' parameter must be numeric between 0 and 1 ([0,1]).')
        }
        if (any(alphas > 0.05)) {
            warning(
                'Most journals and scientists consider an alpha larger than 0.05 untrustworthy. We therefore recommend using alpha values smaller than or queal to 0.05'
            )
        }
        if (!featureToExtract %in%
            c('isoformUsage', 'isoformExp', 'geneExp', 'all')
        ) {
            stop(
                'The \'featureToExtract\' argument must be  \'isoformUsage\', \'isoformExp\',  \'geneExp\' or \'all\''
            )
        }

        ### Check annotation
        okAnnot <-
            c(
                'ORF',
                'NMD_status',
                'coding_potential',
                'signal_peptide_identified',
                'domains_identified',
                'intron_retention',
                'switch_consequences',
                'isoform_class_code'
            )
        if ('all' %in% annotationToAnalyze) {
            annotationToAnalyze <- okAnnot
        }
        if (!all(annotationToAnalyze %in% okAnnot)) {
            stop(
                paste(
                    'The \'annotationToAnalyze\' argument must be a one (or multiple) of: \'',
                    paste(annotationToAnalyze, collapse = '\', \''),
                    '\'',
                    sep = ''
                )
            )
        }

        ### Replace with annotaion names with collum names
        annotToExtract <-
            unique(gsub('NMD_status|ORF', 'PTC', annotationToAnalyze))
        annotToExtract <- gsub(
            'coding_potential',
            'codingPotential',
            annotToExtract
        )
        annotToExtract <- gsub(
            'domains_identified'  ,
            'domain_identified',
            annotToExtract
        )
        annotToExtract <-gsub(
            'intron_retention',
            'IR',
            annotToExtract
        )
        annotToExtract <- gsub(
            'switch_consequences' ,
            'switchConsequencesGene',
            annotToExtract
        )
        annotToExtract <- gsub(
            'isoform_class_code',
            'class_code',
            annotToExtract
        )

        if (length(annotToExtract) == 0) {
            stop(
                'Somthing in the annoation decoding went wrong - please send a small example dataset reconstructing the mistake to the developers.'
            )
        }
    }

    ### Extract annotation
    if (TRUE) {
        switchAnalyzeRlist$isoformFeatures$comparison <- paste(
            switchAnalyzeRlist$isoformFeatures$condition_1,
            switchAnalyzeRlist$isoformFeatures$condition_2,
            sep = ' vs '
        )

        isoformsToAnalyze <- extractSigData(
            switchAnalyzeRlist = switchAnalyzeRlist,
            alpha = alpha,
            dIFcutoff = dIFcutoff,
            log2FCcutoff = log2FCcutoff,
            featureToExtract = featureToExtract
        )

        colToExtract <-
            c('iso_ref',
              'isoform_id',
              'comparison',
              'IF1',
              'IF2',
              annotToExtract)
        colToExtract <-
            intersect(colToExtract,
                      colnames(switchAnalyzeRlist$isoformFeatures))

        isoData <- switchAnalyzeRlist$isoformFeatures[
            which(
                switchAnalyzeRlist$isoformFeatures$iso_ref %in%
                    isoformsToAnalyze
            ),
            na.omit(match(
                colToExtract ,
                colnames(switchAnalyzeRlist$isoformFeatures)
            ))
            ]



        if (nrow(isoData) == 0) {
            stop('No data left after filtering')
        }
    }

    ### Overwrite annotation with propper categories
    if (TRUE) {
        # ORF
        if (!is.null(isoData$PTC) &
            'ORF' %in% annotationToAnalyze) {
            # ORF
            isoData$ORF <- 'With ORF'
            isoData$ORF[which(is.na(isoData$PTC))] <- 'Without ORF'
        }
        # PTC
        if (!is.null(isoData$PTC) &
            'NMD_status' %in% annotationToAnalyze) {
            # PTC
            isoData$PTC <-
                ifelse(test = isoData$PTC,
                       yes = 'NMD sensitive',
                       no = 'NMD insensitive')
        } else {
            isoData$PTC <- NULL
        }
        # coding potential
        if (!is.null(isoData$codingPotential)) {
            isoData$codingPotential <-
                ifelse(test = isoData$codingPotential,
                       yes = 'Isoform is coding',
                       no = 'Isoform is non-coding')
        }
        # signal peptide
        if (!is.null(isoData$signal_peptide_identified)) {
            isoData$signal_peptide_identified <-
                ifelse(
                    test = isoData$signal_peptide_identified == 'yes',
                    yes = 'With signal peptide',
                    no = 'Without signal peptide'
                )
        }
        # protein domains
        if (!is.null(isoData$domain_identified)) {
            isoData$domain_identified <-
                ifelse(
                    test = isoData$domain_identified == 'yes',
                    yes = 'With protein domain',
                    no = 'Without protein domain'
                )
        }
        # intron retention
        if (!is.null(isoData$IR)) {
            isoData$IR <-
                ifelse(test = isoData$IR > 0,
                       yes = 'With intron retention',
                       no = 'Without intron retention')
        }
        # switch consequences
        if (!is.null(isoData$switchConsequences)) {
            isoData$switchConsequences <-
                ifelse(
                    test = isoData$switchConsequences,
                    yes = 'With switch consequence',
                    no = 'Without switch consequence'
                )
        }
        # class code
        if (!is.null(isoData$class_code)) {
            isoData$class_code <-
                paste('class code: \"', isoData$class_code, '\"', sep = '')
        }
    }

    ### Prepare for plotting
    if (TRUE) {
        # melt categories
        isoDataMelt <-
            reshape2::melt(isoData,
                 id.vars = c(
                     'iso_ref', 'isoform_id', 'comparison', 'IF1', 'IF2'
                 ))
        isoDataMelt <-
            isoDataMelt[which(!is.na(isoDataMelt$value)), ]

        # melt IF
        colnames(isoDataMelt)[match(
            c('variable', 'value'), colnames(isoDataMelt)
        )] <- c('category', 'isoform_feature')
        isoDataMelt <-
            reshape2::melt(
                isoDataMelt,
                id.vars = c(
                    'iso_ref',
                    'isoform_id',
                    'comparison',
                    'category',
                    'isoform_feature'
                )
            )

        # massage category
        isoDataMelt$category <- as.character(isoDataMelt$category)

        isoDataMelt$category <- gsub(
            'PTC',
            'NMD Status',
            isoDataMelt$category
        )
        isoDataMelt$category <-
            gsub('codingPotential',
                 'Coding Potential',
                 isoDataMelt$category)
        isoDataMelt$category <-
            gsub('signal_peptide_identified',
                 'Signal Peptide',
                 isoDataMelt$category)
        isoDataMelt$category <-
            gsub('domain_identified',
                 'Protein Domains',
                 isoDataMelt$category)
        isoDataMelt$category <-
            gsub('IR','Intron Retention',isoDataMelt$category)
        isoDataMelt$category <-
            gsub('switchConsequences',
                 'Switch Consequence',
                 isoDataMelt$category)
        isoDataMelt$category <-
            gsub('class_code','Isoform Class',isoDataMelt$category)

        isoDataMelt$comparison2 <-
            paste(isoDataMelt$comparison, '\n(IF1 vs IF2)', sep = '')
    }

    ### Calculate statistics
    if (TRUE) {
        mySigTest <-
            plyr::ddply(
                isoDataMelt,
                .variables = c('comparison', 'category', 'isoform_feature'),
                .drop = TRUE,
                .fun = function(aDF) {
                    data1 <- aDF$value[which(aDF$variable == 'IF1')]
                    data2 <- aDF$value[which(aDF$variable == 'IF2')]

                    myTest <- suppressWarnings(wilcox.test(data1,
                                                           data2))

                    myResult <- data.frame(
                        n = length(data1),
                        medianIF1 = median(data1),
                        medianIF2 = median(data2)
                    )
                    myResult$medianDIF <-
                        myResult$medianIF2 - myResult$medianIF1

                    myResult$wilcoxPval <- myTest$p.value
                    myResult$ymax <- max(c(data1, data2))

                    return(myResult)
                }
            )
        mySigTest$ymax <- mySigTest$ymax * 1.05
        mySigTest$wilcoxQval <-
            p.adjust(mySigTest$wilcoxPval, method = 'fdr')
        mySigTest$significance <-
            sapply(mySigTest$wilcoxQval, function(x)
                evalSig(x, alphas))


        mySigTest <-
            plyr::ddply(
                mySigTest,
                .variables = 'category',
                .fun = function(aDF) {
                    aDF$isoform_feature <- factor(aDF$isoform_feature)
                    aDF$idNr <- as.numeric(aDF$isoform_feature)
                    return(aDF)
                }
            )

        mySigTest$comparison2 <-
            paste(mySigTest$comparison, '\n(IF1 vs IF2)', sep = '')

    }

    ### Plot result
    if (plot) {
        # start plot
        if (violinPlot) {
            p1 <- ggplot() +
                geom_violin(
                    data = isoDataMelt,
                    aes(
                        x = isoform_feature,
                        y = value,
                        fill = variable
                    ),
                    scale = 'area'
                ) +
                stat_summary(
                    data = isoDataMelt,
                    aes(
                        x = isoform_feature,
                        y = value,
                        fill = variable
                    ),
                    fun.y = medianQuartile,
                    geom = 'point',
                    position = position_dodge(width = 0.9),
                    size = 2
                )
        } else {
            p1 <- ggplot() +
                geom_boxplot(data = isoDataMelt, aes(
                    x = isoform_feature,
                    y = value,
                    fill = variable
                ))
        }

        # add significance
        p1 <- p1 +
            geom_text(
                data = mySigTest,
                aes(x = isoform_feature, y = ymax, label = significance),
                vjust = -0.2,
                size = localTheme$text$size * 0.3
            ) +
            geom_segment(data = mySigTest, aes(
                x = idNr - 0.25,
                xend = idNr + 0.25,
                y = ymax,
                yend = ymax
            ))

        # build rest of plot
        p1 <- p1 +
            facet_grid(comparison2 ~ category,
                       scales = 'free_x',
                       space = 'free_x') +
            localTheme +
            theme(strip.text.y = element_text(angle = 0)) +
            theme(axis.text.x = element_text(
                angle = -45,
                hjust = 0,
                vjust = 1
            )) +
            scale_fill_discrete(name = NULL) + theme(legend.position = "top") +
            labs(x = 'Isoform feature', y = 'Isoform Usage (IF)') +
            #coord_cartesian(ylim=c(0,1.25))
            coord_cartesian(ylim = c(0, 1.1 + max(c(
                0, 0.02 * (length(unique(
                    mySigTest$comparison2
                )) - 2)
            ))))

        print(p1)
    }

    ### Return result
    if (returnResult) {
        mySigTest2 <-
            mySigTest[, c(
                'comparison',
                'category',
                'isoform_feature',
                'n',
                'medianIF1',
                'medianIF2',
                'medianDIF',
                'wilcoxPval',
                'wilcoxQval',
                'significance'
            )]
        mySigTest2 <-
            mySigTest2[order(mySigTest2$comparison,
                             mySigTest2$category,
                             mySigTest2$isoform_feature), ]
        return(mySigTest2)
    }
}

# for backward compatability
extractGenomeWideAnalysis <- function(
    switchAnalyzeRlist,
    featureToExtract = 'isoformUsage',
    annotationToAnalyze = 'all',
    alpha=0.05,
    dIFcutoff = 0.1,
    log2FCcutoff = 1,
    violinPlot=TRUE,
    alphas=c(0.05, 0.001),
    localTheme=theme_bw(),
    plot=TRUE,
    returnResult=TRUE
) {
    extractConsequenceGenomeWide(
        switchAnalyzeRlist=switchAnalyzeRlist,
        featureToExtract=featureToExtract,
        annotationToAnalyze=annotationToAnalyze,
        alpha=alpha,
        dIFcutoff=dIFcutoff,
        log2FCcutoff=log2FCcutoff,
        violinPlot=violinPlot,
        alphas=alphas,
        localTheme=localTheme,
        plot=plot,
        returnResult=returnResult
    )
}
