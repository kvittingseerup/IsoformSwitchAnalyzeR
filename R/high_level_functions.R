isoformSwitchAnalysisPart1 <- function(
    switchAnalyzeRlist,
    alpha = 0.05,
    dIFcutoff = 0.1,
    switchTestMethod = 'DEXSeq',
    orfMethod = 'longest',
    genomeObject = NULL,
    cds = NULL,
    pathToOutput = getwd(),
    outputSequences = TRUE,
    prepareForWebServers = FALSE, # to keep backcompatability
    overwriteORF = FALSE,
    quiet = FALSE
) {
    isConditional <- switchAnalyzeRlist$sourceId != 'preDefinedSwitches'
    hasQuant <- ! all(is.na(switchAnalyzeRlist$isoformFeatures$IF_overall))

    nrAnalysis <- 3

    ### Identfy input type
    if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist')        {
        stop(
            'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
        )
    }

    ### Test input
    if(TRUE) {
        ntAlreadyInSwitchList <- ! is.null(switchAnalyzeRlist$ntSequence)
        if( ! ntAlreadyInSwitchList ) {
            if( is.null(genomeObject)) {
                stop(paste(
                    'Since the switchAnalyzeRlist does not contain the',
                    'transcript sequences the BSgenome argument must be used.',
                    '\nIf you used importRdata() to generate the switchAnalyzeRlist',
                    'consider using the isoformNtFasta argument instead (recomended).',
                    sep=' '
                ))
            } else if (class(genomeObject) != 'BSgenome') {
                stop('The genomeObject argument must be a BSgenome')
            }
        }
    }


    ### Run preFilter
    if(hasQuant) {
        switchAnalyzeRlist <-
            preFilter(
                switchAnalyzeRlist = switchAnalyzeRlist,
                removeSingleIsoformGenes = TRUE,
                quiet = TRUE
            )
    }

    ### Test isoform switches
    if(TRUE) {
        if(isConditional) {
            if (switchTestMethod == 'DEXSeq') {
                if (!quiet) {
                    message(
                        paste(
                            'Step 1 of',
                            nrAnalysis,
                            ': Detecting isoform switches...',
                            sep = ' '
                        )
                    )
                }
                switchAnalyzeRlist <-
                    isoformSwitchTestDEXSeq(
                        switchAnalyzeRlist,
                        reduceToSwitchingGenes = TRUE,
                        alpha = alpha,
                        dIFcutoff = dIFcutoff,
                        quiet = TRUE
                    )

            } else if (switchTestMethod == 'DRIMSeq') {
                if (!quiet) {
                    message(
                        paste(
                            'Step 1 of',
                            nrAnalysis,
                            ': Detecting isoform switches (this may take a while)...',
                            sep = ' '
                        )
                    )
                }
                switchAnalyzeRlist <-
                    isoformSwitchTestDRIMSeq(
                        switchAnalyzeRlist,
                        reduceToSwitchingGenes = TRUE,
                        alpha = alpha,
                        dIFcutoff = dIFcutoff,
                        quiet = TRUE
                    )

            } else if (switchTestMethod != 'none') {
                if (!quiet) {
                    message(
                        paste(
                            'Step 1 of',
                            nrAnalysis,
                            ': No isoform switch detection was performed...',
                            sep = ' '
                        )
                    )
                }
                tmp <- 'doNothing'
            } else {
                stop(
                    'Something went wrong with the switch selection - please check input for switch choisce.'
                )
            }


            if (nrow(switchAnalyzeRlist$isoformSwitchAnalysis) == 0) {
                stop('No isoform switches were identified with the current cutoffs.')
            }
        } else {
            if (!quiet) {
                message(
                    paste(
                        'Step 1 of',
                        nrAnalysis,
                        ': No isoform switch detection was performed (since isoform switches were pre-defined)...',
                        sep = ' '
                    )
                )
            }
        }

    }

    ### Predict ORF
    if (all(names(switchAnalyzeRlist) != 'orfAnalysis') | overwriteORF) {
        if (!quiet) {
            message(paste(
                'Step 2 of',
                nrAnalysis,
                ': Predicting open reading frames',
                sep = ' '
            ))
        }

        switchAnalyzeRlist <-
            analyzeORF(
                switchAnalyzeRlist = switchAnalyzeRlist,
                genomeObject = genomeObject,
                cds = cds,
                orfMethod = orfMethod,
                quiet = TRUE
            )
    }

    ### Extract and write sequences
    if (!quiet) {
        message(
            paste(
                'Step 3 of',
                nrAnalysis,
                ': Extracting (and outputting) sequences',
                sep = ' '
            )
        )
    }
    switchAnalyzeRlist <- extractSequence(
        switchAnalyzeRlist = switchAnalyzeRlist,
        genomeObject = genomeObject,
        onlySwitchingGenes = TRUE,
        extractNTseq = TRUE,
        extractAAseq = TRUE,
        addToSwitchAnalyzeRlist = TRUE,
        writeToFile = outputSequences,
        removeLongAAseq = prepareForWebServers,
        alsoSplitFastaFile = prepareForWebServers,
        pathToOutput = pathToOutput,
        quiet = TRUE
    )

    ### Print summary
    if (!quiet) {
        message('\nThe number of isoform switches found were:')
        print(
            extractSwitchSummary(
                switchAnalyzeRlist,
                alpha = alpha,
                dIFcutoff = dIFcutoff
            )
        )
        if(outputSequences) {
            message(
                paste(
                    'The nucleotide and amino acid sequences of these isoforms',
                    'have been outputted to the supplied directory.',
                    '\nThese sequences enabling external analysis of',
                    'protein domians (Pfam), coding potential (CPAT/CPC2) or',
                    'signal peptides (SignalIP).',
                    '\nSee ?analyzeCPAT, ?analyzeCPC2, ?analyzePFAM or?analyzeSignalIP',
                    '(under details) for suggested ways of running these three tools.',
                    sep = ' '
                )
            )
        }

    }
    return(switchAnalyzeRlist)
}

isoformSwitchAnalysisPart2 <- function(
    switchAnalyzeRlist,
    alpha = 0.05,
    dIFcutoff = 0.1,
    n = Inf,
    codingCutoff = NULL,
    removeNoncodinORFs,
    pathToCPATresultFile = NULL,
    pathToCPC2resultFile = NULL,
    pathToPFAMresultFile = NULL,
    pathToIUPred2AresultFile = NULL,
    pathToNetSurfP2resultFile = NULL,
    pathToSignalPresultFile = NULL,
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
    pathToOutput = getwd(),
    fileType = 'pdf',
    asFractionTotal = FALSE,
    outputPlots = TRUE,
    quiet = FALSE
) {
    if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist')        {
        stop(
            'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
        )
    }

    ### Test input
    if (TRUE) {
        if (!is.null(pathToCPATresultFile) & !is.null(pathToCPC2resultFile) ) {
            stop(
                paste(
                    'Since CPC2 and CPAT performs the same type of analysis results should only be suppled to ONE of the \'pathToCPATresultFile\' and \'pathToCPC2resultFile\' arguments.',
                    sep = ' '
                )
            )
        }


        if (!is.null(pathToCPATresultFile)) {
            if (is.null(codingCutoff)) {
                stop(
                    paste(
                        'A cutoff must be supplied to codingCutoff if a CPAT analysis are added to pathToCPATresultFile.',
                        'The cutoff is dependent on spieces analyzed',
                        'see ?analyzeCPAT for more information.',
                        sep = ' '
                    )
                )
            }
        }

        if (!is.null(pathToCPC2resultFile)) {
            if (is.null(codingCutoff)) {
                codingCutoff <- 0.5
            }
        }

        if (!is.null(pathToNetSurfP2resultFile) & !is.null(pathToIUPred2AresultFile) ) {
            stop(
                paste(
                    'Since NetSurfP2 and IUPred2A performs the same type of analysis results should only be suppled to ONE of the \'pathToIUPred2AresultFile\' and \'pathToNetSurfP2resultFile\' arguments.',
                    sep = ' '
                )
            )
        }

    }

    nrAnalysis <-
        2 + as.integer(
            any(c('all', 'intron_retention') %in% consequencesToAnalyze)
        ) + 2 * as.integer(outputPlots)
    analysisDone <- 1

    ### Add annoation
    if (!quiet) {
        message(
            paste(
                'Step',
                analysisDone,
                'of',
                nrAnalysis,
                ': Importing external sequence analysis...',
                sep = ' '
            )
        )
    }

    if (!is.null(pathToCPATresultFile)) {
        switchAnalyzeRlist <-
            analyzeCPAT(
                switchAnalyzeRlist = switchAnalyzeRlist,
                pathToCPATresultFile = pathToCPATresultFile,
                codingCutoff = codingCutoff,
                removeNoncodinORFs = removeNoncodinORFs,
                quiet = TRUE
            )
    }
    if (!is.null(pathToCPC2resultFile)) {
        switchAnalyzeRlist <-
            analyzeCPC2(
                switchAnalyzeRlist = switchAnalyzeRlist,
                pathToCPC2resultFile = pathToCPC2resultFile,
                codingCutoff = codingCutoff,
                removeNoncodinORFs = removeNoncodinORFs,
                quiet = TRUE
            )
    }
    if (!is.null(pathToPFAMresultFile)) {
        switchAnalyzeRlist <-
            analyzePFAM(
                switchAnalyzeRlist = switchAnalyzeRlist,
                pathToPFAMresultFile = pathToPFAMresultFile,
                quiet = TRUE
            )
    }
    if (!is.null(pathToIUPred2AresultFile)) {
        switchAnalyzeRlist <-
            analyzeIUPred2A(
                switchAnalyzeRlist = switchAnalyzeRlist,
                pathToIUPred2AresultFile = pathToIUPred2AresultFile,
                quiet = TRUE
            )
    }
    if (!is.null(pathToNetSurfP2resultFile)) {
        switchAnalyzeRlist <-
            analyzeNetSurfP2(
                switchAnalyzeRlist = switchAnalyzeRlist,
                pathToNetSurfP2resultFile = pathToNetSurfP2resultFile,
                quiet = TRUE
            )
    }
    if (!is.null(pathToSignalPresultFile)) {
        switchAnalyzeRlist <-
            analyzeSignalP(
                switchAnalyzeRlist = switchAnalyzeRlist,
                pathToSignalPresultFile = pathToSignalPresultFile,
                quiet = TRUE
            )
    }
    analysisDone <- analysisDone + 1

    ### Predict intron retentions
    if (any(c('all', 'intron_retention') %in% consequencesToAnalyze)) {
        if (!quiet) {
            message(
                paste(
                    'Step',
                    analysisDone,
                    'of',
                    nrAnalysis,
                    ': Analyzing alternative splicing...',
                    sep = ' '
                )
            )
        }

        switchAnalyzeRlist <-
            analyzeAlternativeSplicing(
                switchAnalyzeRlist = switchAnalyzeRlist,
                onlySwitchingGenes = TRUE,
                alpha = alpha,
                dIFcutoff = dIFcutoff,
                quiet = TRUE
            )
        analysisDone <- analysisDone + 1
    }

    ### Predict functional consequences
    if (!quiet) {
        message(
            paste(
                'Step',
                analysisDone,
                'of',
                nrAnalysis,
                ': Prediciting functional consequences...',
                sep = ' '
            )
        )
    }

    switchAnalyzeRlist <-
        analyzeSwitchConsequences(
            switchAnalyzeRlist = switchAnalyzeRlist,
            consequencesToAnalyze = consequencesToAnalyze,
            alpha = alpha,
            dIFcutoff = dIFcutoff,
            quiet = TRUE
        )
    analysisDone <- analysisDone + 1

    ### Make isoform switch plots
    if (outputPlots) {
        if (!quiet) {
            message(
                paste(
                    'Step',
                    analysisDone,
                    'of',
                    nrAnalysis,
                    ': Making indidual isoform switch plots...',
                    sep = ' '
                )
            )
        }

        switchPlotTopSwitches(
            switchAnalyzeRlist = switchAnalyzeRlist,
            alpha = alpha,
            dIFcutoff = dIFcutoff,
            n = n,
            pathToOutput = pathToOutput,
            filterForConsequences = TRUE,
            splitFunctionalConsequences = FALSE,
            fileType = fileType,
            quiet = TRUE
        )

        analysisDone <- analysisDone + 1
    }

    ### Make overall consequences
    if (outputPlots) {
        if (!quiet) {
            message(
                paste(
                    'Step',
                    analysisDone,
                    'of',
                    nrAnalysis,
                    ': Analyzing combined consequences plot...',
                    sep = ' '
                )
            )
        }

        ### Summary
        if(TRUE) {
            if (fileType == 'pdf') {
                pdf(
                    file = paste(
                        pathToOutput,
                        '/common_switch_consequences.pdf',
                        sep = ''
                    ),
                    width = 10,
                    height = 7
                )
            } else {
                png(
                    filename = paste(
                        pathToOutput,
                        '/common_switch_consequences.png',
                        sep = ''
                    ),
                    width = 1000,
                    height = 700
                )
            }
            extractConsequenceSummary(
                switchAnalyzeRlist = switchAnalyzeRlist,
                asFractionTotal = asFractionTotal,
                alpha = alpha,
                dIFcutoff = dIFcutoff,
                plot = TRUE,
                returnResult = FALSE
            )
            dev.off()

        }

        ### Consequence enrichment
        if(TRUE) {
            ### test numbers found
            testNumbers <- extractConsequenceEnrichment(
                switchAnalyzeRlist = switchAnalyzeRlist,
                alpha = alpha,
                dIFcutoff = dIFcutoff,
                plot=FALSE,
                returnResult = TRUE
            )

            enougthEvents <- any(
                testNumbers$nUp + testNumbers$nDown >= 10
            )

            if(enougthEvents) {
                if (fileType == 'pdf') {
                    pdf(
                        file = paste(
                            pathToOutput,
                            '/switch_consequences_enrichment.pdf',
                            sep = ''
                        ),
                        width = 10,
                        height = 7
                    )
                } else {
                    png(
                        filename = paste(
                            pathToOutput,
                            '/switch_consequences_enrichment.png',
                            sep = ''
                        ),
                        width = 1000,
                        height = 700
                    )
                }
                extractConsequenceEnrichment(
                    switchAnalyzeRlist = switchAnalyzeRlist,
                    alpha = alpha,
                    dIFcutoff = dIFcutoff,
                    plot=TRUE,
                    returnResult = FALSE
                )
                dev.off()
            } else {
                warning(
                    'extractConsequenceEnrichment() could not be run because to few consequences were detected.'
                )
            }


        }

        ### Splicing enrichment
        if(TRUE) {
            testNumbers <- extractSplicingEnrichment(
                switchAnalyzeRlist = switchAnalyzeRlist,
                alpha = alpha,
                dIFcutoff = dIFcutoff,
                plot=FALSE,
                returnResult = TRUE
            )

            enougthEvents <- any(
                testNumbers$nUp + testNumbers$nDown >= 10
            )

            if(enougthEvents) {
                if (fileType == 'pdf') {
                    pdf(
                        file = paste(
                            pathToOutput,
                            '/splicing_enrichment.pdf',
                            sep = ''
                        ),
                        width = 10,
                        height = 7
                    )
                } else {
                    png(
                        filename = paste(
                            pathToOutput,
                            '/splicing_enrichment.png',
                            sep = ''
                        ),
                        width = 1000,
                        height = 700
                    )
                }
                extractSplicingEnrichment(
                    switchAnalyzeRlist = switchAnalyzeRlist,
                    alpha = alpha,
                    dIFcutoff = dIFcutoff,
                    plot=TRUE,
                    returnResult = FALSE
                )
                dev.off()
            } else {
                warning(
                    'extractSplicingEnrichment() could not be run because to few splicing differences were detected.'
                )
            }


        }
    }

    ### Print summary
    if (!quiet) {
        message(
            '\nThe number of isoform switches with functional consequences identified were:'
        )
        print(
            extractSwitchSummary(
                switchAnalyzeRlist,
                alpha = alpha,
                dIFcutoff = dIFcutoff,
                filterForConsequences = TRUE
            )
        )
        if (outputPlots) {
            message(
                paste(
                    'The switch analysis plot for each of these, as well as a plot summarizing the functional consequences',
                    'have been outputted to the folder specified by \'pathToOutput\'.',
                    sep = ' '
                )
            )
        }
    }

    return(switchAnalyzeRlist)
}

isoformSwitchAnalysisCombined <- function(
    switchAnalyzeRlist,
    alpha = 0.05,
    dIFcutoff = 0.1,
    switchTestMethod = 'DEXSeq',
    n = NA,
    pathToOutput = getwd(),
    overwriteORF = FALSE,
    outputSequences = FALSE,
    genomeObject,
    orfMethod = 'longest',
    cds = NULL,
    consequencesToAnalyze = c('intron_retention', 'ORF_seq_similarity', 'NMD_status'),
    fileType = 'pdf',
    asFractionTotal = FALSE,
    outputPlots = TRUE,
    quiet = FALSE
) {
    ### Run part 1
    if (!quiet) {
        message('\nPART 1: EXTRACTING ISOFORM SWTICH SEQUENCES')
    }
    switchAnalyzeRlist <-
        isoformSwitchAnalysisPart1(
            switchAnalyzeRlist = switchAnalyzeRlist,
            alpha = alpha,
            dIFcutoff = dIFcutoff,
            switchTestMethod = switchTestMethod,
            pathToOutput = pathToOutput,
            genomeObject = genomeObject,
            orfMethod = orfMethod,
            cds = cds,
            outputSequences = outputSequences,
            overwriteORF = overwriteORF,
            quiet = quiet
        )

    ### Run part 2 without annoation
    if (!quiet) {
        message('\nPART 2: PLOTTING ISOFORM SEQUENCES')
    }
    switchAnalyzeRlist <-
        isoformSwitchAnalysisPart2(
            switchAnalyzeRlist = switchAnalyzeRlist,
            alpha = alpha,
            dIFcutoff = dIFcutoff,
            n = n,
            consequencesToAnalyze = consequencesToAnalyze,
            pathToOutput = pathToOutput,
            fileType = fileType,
            asFractionTotal = asFractionTotal,
            outputPlots = outputPlots,
            quiet = quiet
        )

    return(switchAnalyzeRlist)
}
