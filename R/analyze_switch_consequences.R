analyzeSwitchConsequences <- function(
    switchAnalyzeRlist,
    consequencesToAnalyze = c(
        'intron_retention',
        'coding_potential',
        'ORF_seq_similarity',
        'NMD_status',
        'domains_identified',
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
            'signal_peptide_identified'
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
            if (is.null(switchAnalyzeRlist$intronRetentionAnalysis)) {
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

    }

    if (showProgress & !quiet) {
        progressBar <- 'text'
    } else {
        progressBar <- 'none'
    }

    ### Subset to relevant data
    if (!quiet) {
        message('Step 1 of 4: Extracting genes with isoform switches...')
    }
    if (TRUE) {
        localData <- switchAnalyzeRlist$isoformFeatures[which(
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value < alpha &
                abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
        ),
        c(
            'iso_ref',
            'gene_ref',
            'isoform_switch_q_value',
            'gene_switch_q_value',
            'dIF'
        )]
        if (!nrow(localData)) {
            stop('No genes were considered switching with the used cutoff values')
        }

        ### add switch direction
        localData$switchDirection <- NA
        localData$switchDirection[which(sign(localData$dIF) ==  1)] <- 'up'
        localData$switchDirection[which(sign(localData$dIF) == -1)] <- 'down'

        ### split based on genes and conditions
        localDataList <-
            split(localData, f = localData$gene_ref, drop = TRUE)

        ### Extract pairs of isoforms passing the filters
        pairwiseIsoComparisonList <-
            llply(
                .data = localDataList,
                .progress = 'none',
                .fun = function(aDF) {
                    # aDF <- localDataList[[171]]
                    isoResTest <-
                        any(!is.na(
                            switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
                        ))
                    if (isoResTest) {
                        sigIso <- aDF$iso_ref[which(
                            aDF$isoform_switch_q_value < alpha &
                                abs(aDF$dIF) > dIFcutoff
                        )]
                    } else {
                        sigIso <- aDF$iso_ref[which(
                            aDF$gene_switch_q_value < alpha &
                                abs(aDF$dIF) > dIFcutoff
                        )]
                    }
                    if (length(sigIso) == 0) {
                        return(NULL)
                    }

                    ### reduce to significant if nessesary
                    if (onlySigIsoforms) {
                        aDF <- aDF[which(aDF$iso_ref %in% sigIso), ]
                    }
                    if (nrow(aDF) < 2) {
                        return(NULL)
                    }

                    ### make sure there are both up and down
                    if (!all(c('up', 'down') %in% aDF$switchDirection)) {
                        return(NULL)
                    }

                    ### extract pairs of isoforms
                    upIso   <-
                        as.vector(aDF$iso_ref[which(
                            aDF$switchDirection == 'up'
                        )])
                    downIso <-
                        as.vector(aDF$iso_ref[which(
                            aDF$switchDirection == 'down'
                        )])

                    allIsoCombinations <-
                        setNames(
                            base::expand.grid(
                                upIso,
                                downIso,
                                stringsAsFactors = FALSE,
                                KEEP.OUT.ATTRS = FALSE
                            ),
                            nm = c('iso_ref_up', 'iso_ref_down')
                        )

                    ### Reduce to those where at least one of them is significant
                    allIsoCombinations <-
                        allIsoCombinations[which(
                            allIsoCombinations$iso_ref_up %in% sigIso |
                                allIsoCombinations$iso_ref_down %in% sigIso
                        ), ]

                    ### Add gen ref
                    allIsoCombinations$gene_ref    <- aDF$gene_ref[1]

                    return(allIsoCombinations)
                }
            )

        ### Remove empty entries
        pairwiseIsoComparisonList <-
            pairwiseIsoComparisonList[which(
                ! sapply(pairwiseIsoComparisonList, is.null)
            )]
        if (length(pairwiseIsoComparisonList) == 0) {
            stop('No candidate genes with the required cutoffs were found')
        }

        ### Conver to data.frame
        pairwiseIsoComparison <-
            myListToDf(pairwiseIsoComparisonList, addOrignAsColumn = FALSE)

        ### Add isoform names
        pairwiseIsoComparison$isoformUpregulated   <-
            switchAnalyzeRlist$isoformFeatures$isoform_id[match(
                pairwiseIsoComparison$iso_ref_up,
                switchAnalyzeRlist$isoformFeatures$iso_ref
            )]
        pairwiseIsoComparison$isoformDownregulated <-
            switchAnalyzeRlist$isoformFeatures$isoform_id[match(
                pairwiseIsoComparison$iso_ref_down,
                switchAnalyzeRlist$isoformFeatures$iso_ref
            )]

        ### Extract paris to compare
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
    if (TRUE) {
        consequencesOfIsoformSwitching <- dlply(
            .data = pairwiseIsoComparisonUniq,
            .variables = 'comparison',
            .parallel = FALSE,
            .inform = TRUE,
            .progress = progressBar,
            .fun = function(aDF) {
                # aDF <- pairwiseIsoComparisonUniq[1,]
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
    if (!quiet) {
        message(paste(
            'Step 3 of 4: Massaging isoforms comparisons results...',
            sep = ' '
        ))
    }
    if (TRUE) {
        ### Convert from list to df
        consequencesOfIsoformSwitchingDf <-
            myListToDf(consequencesOfIsoformSwitching)

        ### Add the origin info
        consequencesOfIsoformSwitchingDfcomplete <- merge(
            consequencesOfIsoformSwitchingDf,
            pairwiseIsoComparison,
            by = c('isoformUpregulated', 'isoformDownregulated')
        )
        ### Add additional information
        consequencesOfIsoformSwitchingDfcomplete <- merge(
            consequencesOfIsoformSwitchingDfcomplete,
            switchAnalyzeRlist$isoformFeatures[match(
                unique(consequencesOfIsoformSwitchingDfcomplete$gene_ref),
                switchAnalyzeRlist$isoformFeatures$gene_ref
            ),
            c('gene_ref',
              'gene_id',
              'gene_name',
              'condition_1',
              'condition_2')],
            by = 'gene_ref'
        )

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
    if (!quiet) {
        message('Step 4 of 4: Preparing output...')
    }
    if (TRUE) {
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
            'signal_peptide_identified'
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
            if (is.null(switchAnalyzeRlist$intronRetentionAnalysis)) {
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
                'signal_peptide_identified'
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
                merge(transcriptData,
                      orfData,
                      by = 'isoform_id',
                      all.x = TRUE)
        }

        ### intron retention data
        if ('intron_retention' %in% consequencesToAnalyze) {
            localIRdata <-
                switchAnalyzeRlist$intronRetentionAnalysis[which(
                    switchAnalyzeRlist$intronRetentionAnalysis$isoform_id %in%
                        isoformsToAnalyze
                ), ]
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
        ### Makre result data.frame
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
                pairwiseAlignment(pattern = upNtSeq,
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

        ### Intron retention by analyzing spliceR annotation
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
                        pairwiseAlignment(
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
                        pairwiseAlignment(
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
                        pairwiseAlignment(
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
                commonDomains <-
                    table(unlist(lapply(domanDataSplit, function(x)
                        x$hmm_name)))
                commonDomains <-
                    commonDomains[which(commonDomains == 2)]

                # Any common domains
                if (length(commonDomains)) {
                    # for ecah domain analyze legth differences
                    domainTruncationStatus <-
                        unique(sapply(names(commonDomains), function(aDomain) {
                            localDomian <-
                                ldply(
                                    domanDataSplit,
                                    .fun = function(x)
                                        x[which(x$hmm_name == aDomain), ]
                                )
                            localDomian$length <-
                                localDomian$orf_aa_end -
                                localDomian$orf_aa_start + 1

                            domainLengthDifferent <-
                                abs(diff(localDomian$length)) > AaCutoff

                            if (!is.null(AaFracCutoff)) {
                                maxIndex <- which.max(localDomian$length)
                                fractionDifference <-
                                    abs(diff(localDomian$length)) /
                                    localDomian$length[maxIndex] > AaFracCutoff
                                domainLengthDifferent <-
                                    domainLengthDifferent & fractionDifference
                            }

                            if (domainLengthDifferent) {
                                if (addDescription) {
                                    shortIso <- localDomian$.id[which.min(
                                        localDomian$length
                                    )]

                                    if (shortIso == upIso) {
                                        domianLegnthStatus <-
                                            'Domain length gain'
                                    } else {
                                        domianLegnthStatus <-
                                            'Domain length loss'
                                    }
                                } else {
                                    domianLegnthStatus <- TRUE
                                }

                            } else {
                                domianLegnthStatus <- NA
                            }

                            return(domianLegnthStatus)
                        }))

                    ### Summarize differences
                    localIndex <-
                        which(isoComparison$featureCompared == 'domain_length')
                    if (all(is.na(domainTruncationStatus))) {
                        isoComparison$isoformsDifferent[localIndex] <- FALSE
                    } else if (length(na.omit(domainTruncationStatus)) == 2) {
                        isoComparison$isoformsDifferent[localIndex] <- TRUE

                        if (addDescription) {
                            isoComparison$switchConsequence[localIndex] <-
                                'Domain length gain and loss'
                        }
                    } else if (length(na.omit(domainTruncationStatus)) == 1) {
                        isoComparison$isoformsDifferent[localIndex] <- TRUE
                        if (addDescription) {
                            isoComparison$switchConsequence[localIndex] <-
                                na.omit(domainTruncationStatus)
                        }
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

    }

    if (onlyRepportDifferent) {
        isoComparison <-
            isoComparison[which(isoComparison$isoformsDifferent), ]
    }

    return(isoComparison)
}

extractConsequenceSummary <- function(
    switchAnalyzeRlist,
    consequencesToAnalyze = 'all',
    includeCombined = FALSE,
    asFractionTotal = FALSE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    plot = TRUE,
    plotGenes = FALSE,
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
        if (asFractionTotal) {
            if (alpha < 0 |
                alpha > 1) {
                warning('The alpha parameter should usually be between 0 and 1 ([0,1]).')
            }
            if (alpha > 0.05) {
                warning(
                    'Most journals and scientists consider an alpha larger than 0.05 untrustworthy. We therefore recommend using alpha values smaller than or queal to 0.05'
                )
            }
        }

        ### Consequences to analyze
        acceptedTypes <- c(
            # Transcript
            'tss',
            'tts',
            'last_exon',
            'isoform_length',
            'isoform_seq_similarity',
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
            'signal_peptide_identified'
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


    }

    startCapitalLetter <- function(aVec) {
        paste(
            toupper( substr(aVec, 1, 1) ),
            substr(aVec, 2, nchar(aVec)),
            sep = ""
        )
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
        ddply(
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
