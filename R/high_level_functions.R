isoformSwitchAnalysisPart1 <- function(
    ### Core arguments
    switchAnalyzeRlist,

    ### Annotation arguments
    genomeObject = NULL,
    pathToGTF = NULL,

    ### Output arguments
    prepareForWebServers,
    pathToOutput = getwd(),
    outputSequences = TRUE,

    ### Other arguments
    quiet = FALSE
) {
    isConditional <- switchAnalyzeRlist$sourceId != 'preDefinedSwitches'

    nrAnalysis <- 3

    ### Identfy input type
    if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist')        {
        stop(
            'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
        )
    }

    ### Test input
    if(TRUE) {
        ### Test sequence annotation
        if( is.null(switchAnalyzeRlist$ntSequence) ) {
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
        } else {
            if( ! is.null(genomeObject)) {
                warning(paste0(
                    'Since isoform sequences are already annotated in the switchAnalyzeRlist',
                    '\n   these will be used (meaning the \'genomeObject\' argument will be ignored.'
                ))
                genomeObject <- NULL
            }
        }

        ### Test ORF annotation
        if( is.null(switchAnalyzeRlist$orfAnalysis) ) {
            if( is.null(pathToGTF)) {
                stop('When ORFs are not annotated you must supply the GTF with known annotation to the \'pathToGTF\' argument.')
            }
        }
    }


    ### preFilter
    if(isConditional) {
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

            if(any( switchAnalyzeRlist$conditions$nrReplicates > 5)) {
                switchAnalyzeRlist <-
                    isoformSwitchTestSatuRn(
                        switchAnalyzeRlist,
                        reduceToSwitchingGenes = TRUE,
                        quiet = TRUE
                    )
            } else {
                switchAnalyzeRlist <-
                    isoformSwitchTestDEXSeq(
                        switchAnalyzeRlist,
                        reduceToSwitchingGenes = TRUE,
                        quiet = TRUE
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
    if ( is.null(switchAnalyzeRlist$orfAnalysis) ) {
        if (!quiet) {
            message(paste(
                'Step 2 of',
                nrAnalysis,
                ': Adding and predicting open reading frames',
                sep = ' '
            ))
        }

        ### Add known annoation
        switchAnalyzeRlist <- addORFfromGTF(
            switchAnalyzeRlist = switchAnalyzeRlist,
            pathToGTF = pathToGTF,
            quiet = TRUE
        )

        ### Predict novel once (if any are missing)
        if ( any( switchAnalyzeRlist$orfAnalysis$orf_origin == 'not_annotated_yet' )) {
            switchAnalyzeRlist <- analyzeNovelIsoformORF(
                switchAnalyzeRlist = switchAnalyzeRlist,
                analysisAllIsoformsWithoutORF = TRUE,
                genomeObject = genomeObject
            )
        }

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
                switchAnalyzeRlist
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
    ### Core arguments
    switchAnalyzeRlist,

    ### External annotation arguments
    codingCutoff = NULL,
    removeNoncodinORFs,
    pathToCPATresultFile = NULL,
    pathToCPC2resultFile = NULL,
    pathToPFAMresultFile = NULL,
    pathToIUPred2AresultFile = NULL,
    pathToNetSurfP2resultFile = NULL,
    pathToSignalPresultFile = NULL,
    pathToDeepLoc2resultFile = NULL,
    pathToDeepTMHMMresultFile = NULL,

    ### Analysis and output arguments
    n = Inf,
    consequencesToAnalyze = c(
        'intron_retention',
        'coding_potential',
        'ORF_seq_similarity',
        'NMD_status',
        'domains_identified',
        'domain_isotype',
        'IDR_identified',
        'IDR_type',
        'signal_peptide_identified'
    ),
    pathToOutput = getwd(),
    fileType = 'pdf',
    outputPlots = TRUE,

    ### Other arguments
    quiet = FALSE
) {
    ### Test input
    if (TRUE) {
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist')        {
            stop(
                'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
            )
        }

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
    if (!is.null(pathToSignalPresultFile)) {
        switchAnalyzeRlist <-
            analyzeDeepLoc2(
                switchAnalyzeRlist = switchAnalyzeRlist,
                pathToDeepLoc2resultFile = pathToDeepLoc2resultFile,
                quiet = TRUE
            )
    }
    if (!is.null(pathToSignalPresultFile)) {
        switchAnalyzeRlist <-
            analyzeDeepTMHMM(
                switchAnalyzeRlist = switchAnalyzeRlist,
                pathToDeepTMHMMresultFile = pathToDeepTMHMMresultFile,
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
                ': Predicting functional consequences...',
                sep = ' '
            )
        )
    }

    switchAnalyzeRlist <-
        analyzeSwitchConsequences(
            switchAnalyzeRlist = switchAnalyzeRlist,
            consequencesToAnalyze = consequencesToAnalyze,
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
            print(
                extractConsequenceSummary(
                    switchAnalyzeRlist = switchAnalyzeRlist,
                    plot = TRUE,
                    returnResult = FALSE
                )
            )
            dev.off()

        }

        ### Consequence enrichment
        if(TRUE) {
            ### test numbers found
            testNumbers <- extractConsequenceEnrichment(
                switchAnalyzeRlist = switchAnalyzeRlist,
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
    ### Core arguments
    switchAnalyzeRlist,

    ### Annotation arguments
    genomeObject = NULL,
    pathToGTF = NULL,

    ### Analysis and output arguments
    n = Inf,
    consequencesToAnalyze = c('intron_retention', 'ORF_seq_similarity', 'NMD_status'),
    pathToOutput = getwd(),
    fileType = 'pdf',
    outputPlots = TRUE,

    ### Other arguments
    quiet = FALSE
) {
    ### Run part 1
    if (!quiet) {
        message('\nPART 1: EXTRACTING ISOFORM SWTICH SEQUENCES')
    }
    switchAnalyzeRlist <-
        isoformSwitchAnalysisPart1(
            switchAnalyzeRlist = switchAnalyzeRlist,
            pathToOutput = pathToOutput,
            genomeObject = genomeObject,
            pathToGTF = pathToGTF,
            outputSequences = FALSE,
            prepareForWebServers = FALSE,
            quiet = quiet
        )

    ### Run part 2 without annoation
    if (!quiet) {
        message('\nPART 2: PLOTTING ISOFORM SEQUENCES')
    }
    switchAnalyzeRlist <-
        isoformSwitchAnalysisPart2(
            switchAnalyzeRlist = switchAnalyzeRlist,
            n = n,
            consequencesToAnalyze = consequencesToAnalyze,
            pathToOutput = pathToOutput,
            fileType = fileType,
            outputPlots = outputPlots,
            quiet = quiet
        )

    return(switchAnalyzeRlist)
}
