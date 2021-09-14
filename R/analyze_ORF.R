########################################################################
########################### switchAnalyzeRlist #########################
######### Analyze coding potential and identify protein domains ########
########################################################################

getCDS <- function(
    selectedGenome,
    repoName
) {
    ### Test input
    if (length(repoName) > 1)
        stop("getCDS: Please supply only one repository")
    if (!any(selectedGenome %in% c("hg19", "hg38", "mm9", "mm10")))
        stop("getCDS: supported genomes are currently: hg19, hg38, mm9, mm10")
    if (!any(repoName %in% c("ensemble", "UCSC", "refseq", "GENCODE")))
        stop("getCDS: Supported repositories are currently: Ensemble (not hg38), UCSC, Refseq, GENCODE (not mm9)")

    ### Make data.frame to translate to tack names
    localRepos <- rbind(
        # hg19
        data.frame(
            genome='hg19',
            repo  =c("ensemble", "UCSC"     , "refseq"         , "GENCODE"),
            track =c("ensGene" , "knownGene", "refGene"        , "wgEncodeGencodeV19"),
            stringsAsFactors = FALSE
        ),
        # hg38
        data.frame(
            genome='hg38',
            repo  =c("ensemble", "UCSC"     , "refseq"         , "GENCODE"),
            track =c(NA        , NA         , 'refSeqComposite', "knownGene"), # gencode is stored as knownGene according to trackNames() annotation
            stringsAsFactors = FALSE
        ),
        # mm9
        data.frame(
            genome='mm9',
            repo  =c("ensemble", "UCSC"     , "refseq"         , "GENCODE"),
            track =c("ensGene" , "knownGene", "refGene"        , NA),
            stringsAsFactors = FALSE
        ),
        # mm10
        data.frame(
            genome='mm10',
            repo  =c("ensemble", "UCSC"     , "refseq"         , "GENCODE"),
            track =c(NA        , "knownGene", "refSeqComposite", "wgEncodeGencode"),
            stringsAsFactors = FALSE
        )
    )
    repoOfInterest <- localRepos[which(
        localRepos$genome == selectedGenome &
            localRepos$repo == repoName
    ),]

    if(is.na(repoOfInterest$track)) {
        stop('Combination of genome and repository is not available')
    }

    ### Make session
    message("Retrieving CDS tables for ", repoOfInterest$repo, "...", sep = "")
    session <- rtracklayer::browserSession("UCSC",url="http://genome-euro.ucsc.edu/cgi-bin/") # KVS : this solves the problem with changing geomes
    GenomeInfoDb::genome(session) <- repoOfInterest$genome
    query <- rtracklayer::ucscTableQuery(session, repoOfInterest$track)

    ### Get data
    cdsTable <- rtracklayer::getTable(query)
    if (repoOfInterest$repo == "ensemble")
        cdsTable <- cdsTable[cdsTable$cdsStartStat != "none",
                             ]
    if (repoOfInterest$repo == "UCSC")
        cdsTable <- cdsTable[cdsTable$cdsStart != cdsTable$cdsEnd,
                             ]
    if (repoOfInterest$repo == "refseq")
        cdsTable <- cdsTable[cdsTable$cdsStart != cdsTable$cdsEnd,
                             ]
    cdsTable <- cdsTable[, c("chrom", "strand", "txStart", "txEnd",
                             "cdsStart", "cdsEnd", "exonCount", "name")]
    message("Retrieved ", nrow(cdsTable), " records...", sep = "")
    utils::flush.console()
    return(new("CDSSet", cdsTable))
}

analyzeORF <- function(
    ### Core arguments
    switchAnalyzeRlist,
    genomeObject = NULL,

    ### Advanced argument
    minORFlength = 100,
    orfMethod = 'longest',
    cds = NULL,
    PTCDistance = 50,
    startCodons = "ATG",
    stopCodons = c("TAA", "TAG", "TGA"),
    showProgress = TRUE,
    quiet = FALSE
) {
    ### check input
    if (TRUE) {
        # Input data
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(paste(
                'The object supplied to \'switchAnalyzeRlist\'',
                'must be a \'switchAnalyzeRlist\''
            ))
        }

        ntAlreadyInSwitchList <- ! is.null(switchAnalyzeRlist$ntSequence)
        if( ! ntAlreadyInSwitchList ) {
            if (class(genomeObject) != 'BSgenome') {
                stop('The genomeObject argument must be a BSgenome')
            }
        }

        # method choice
        if (length(orfMethod) != 1) {
            stop(paste(
                'The \'orfMethod\' must be one of \'mostUpstreamAnnoated\',',
                '\'mostUpstream\', \'longestAnnotated\', \'longest.AnnotatedWhenPossible\' or \'longest\',',
                'not a vector'
            ))
        }
        if (!orfMethod %in% c(
            'mostUpstreamAnnoated',
            'mostUpstream',
            'longest',
            'longestAnnotated',
            'longest.AnnotatedWhenPossible'
            )
        ) {
            stop(paste(
                'The \'orfMethod\' argument must be either of',
                '\'mostUpstreamAnnoated\',\'mostUpstream\',',
                ' \'longestAnnotated\', \'longest.AnnotatedWhenPossible\' or \'longest\' indicating whether',
                'to use the the most upstream annoated start codon or',
                'the longest ORF respectively'
            ))
        }

        if (orfMethod %in% c('mostUpstreamAnnoated', 'longestAnnotated')) {
            if (is.null(cds)) {
                stop(paste(
                    'When using orfMethod is \'mostUpstreamAnnoated\' or',
                    '\'longestAnnotated\', a CDSSet must be supplied to',
                    'the \'cds\' argument '
                ))
            }
        }

    }

    ### Assign paramters
    if (TRUE) {
        useAnnoated <-
            orfMethod %in% c('longestAnnotated', 'mostUpstreamAnnoated','longest.AnnotatedWhenPossible')

        nrStepsToPerform <- 3 + as.numeric(useAnnoated)
        nrStepsPerformed <- 0

        if (showProgress & !quiet) {
            progressBar <- 'text'
        } else {
            progressBar <- 'none'
        }

        myCodonDf <- data.frame(
            codons = c(startCodons, stopCodons),
            meaning = c(
                replicate(n = length(startCodons), expr = 'start'),
                replicate(n = length(stopCodons), expr = 'stop')
            ),
            stringsAsFactors = FALSE
        )

    }

    ### Idenetify overlapping CDS and exons
    if (useAnnoated) {
        if (!quiet) {
            message(
                paste(
                    'Step',
                    nrStepsPerformed + 1,
                    'of',
                    nrStepsToPerform,
                    ': Identifying overlap between supplied CDS and isoforms',
                    sep = ' '
                )
            )
        }

        ### Idenetify overlapping CDS and exons and annoate transcript position
        cds <- unique(cds[, c('chrom', 'strand', 'cdsStart', 'cdsEnd')])

        ### strand specific translation starts
        plusIndex <- which(cds$strand == '+')
        annoatedStartGRangesPlus <-
            unique(GRanges(
                cds$chrom[plusIndex],
                IRanges(
                    start = cds$cdsStart[plusIndex],
                    end = cds$cdsStart[plusIndex]
                ),
                strand = cds$strand[plusIndex]
            ))
        minusIndex <- which(cds$strand == '-')
        annoatedStartGRangesMinus <-
            unique(GRanges(
                cds$chrom[minusIndex],
                IRanges(
                    start = cds$cdsEnd[minusIndex],
                    end = cds$cdsEnd[minusIndex]),
                strand = cds$strand[minusIndex]
            ))

        annoatedStartGRanges <-
            unique(c(annoatedStartGRangesPlus, annoatedStartGRangesMinus))
        annoatedStartGRanges$id <-
            paste('cds_', 1:length(annoatedStartGRanges), sep = '')

        ### Extract exons
        localExons <-  switchAnalyzeRlist$exons[which(
            switchAnalyzeRlist$exons$isoform_id %in%
                switchAnalyzeRlist$isoformFeatures$isoform_id
        ),]
        #localExons <- localExons[which(strand(localExons) %in% c('+','-')),]
        localExons <-
            localExons[which(as.character(localExons@strand) %in% c('+', '-')),]

        localExons <-
            localExons[order(localExons$isoform_id,
                             start(localExons),
                             end(localExons)), ]
        localExons$exon_id <-
            paste('exon_', 1:length(localExons), sep = '')


        ### Find overlaps
        suppressWarnings(
            overlappingAnnotStart <-
                as.data.frame(
                    findOverlaps(
                        query = localExons,
                        subject = annoatedStartGRanges,
                        ignore.strand = FALSE
                    )
                )
        )
        #if (!nrow(overlappingAnnotStart)) {
        #    stop(
        #        'No overlap between CDS and transcripts were found. This is most likely due to a annoation problem around chromosome name.'
        #    )
        #}

        ### Annoate overlap
        overlappingAnnotStart$queryHits <-
            localExons$exon_id[overlappingAnnotStart$queryHits]
        overlappingAnnotStart$subjectHits <-
            annoatedStartGRanges$id[overlappingAnnotStart$subjectHits]
        colnames(overlappingAnnotStart) <- c('exon_id', 'cds_id')
        overlappingAnnotStart$isoform_id <-
            localExons$isoform_id[match(
                overlappingAnnotStart$exon_id, localExons$exon_id
            )]

        ### annoate with genomic start site
        overlappingAnnotStart$cdsGenomicStart <-
            start(annoatedStartGRanges)[match(overlappingAnnotStart$cds_id,
                                              annoatedStartGRanges$id)]


        ## Extract exon information
        myExons <-
            as.data.frame(localExons[which(
                localExons$isoform_id %in% overlappingAnnotStart$isoform_id
            ),])
        myExons <- myExons[sort.list(myExons$isoform_id), ]
        #myExonsSplit <- split(myExons, f=myExons$isoform_id)

        myExonPlus <- myExons[which(myExons$strand == '+'), ]
        myExonPlus$cumSum <-
            unlist(sapply(
                split(myExonPlus$width, myExonPlus$isoform_id),
                function(aVec) {
                    cumsum(c(0, aVec))[1:(length(aVec))]
                }
            ))
        myExonMinus <- myExons[which(myExons$strand == '-'), ]
        myExonMinus$cumSum <-
            unlist(sapply(
                split(myExonMinus$width, myExonMinus$isoform_id),
                function(aVec) {
                    cumsum(c(0, rev(aVec)))[(length(aVec)):1] # reverse
                }
            ))

        myExons2 <- rbind(myExonPlus, myExonMinus)
        myExons2 <-
            myExons2[order(myExons2$isoform_id, myExons2$start, myExons2$end), ]
        #myExonsSplit <- split(myExons2, f=myExons2$isoform_id)

        ### Annoate with exon information
        matchIndex <-
            match(overlappingAnnotStart$exon_id, myExons2$exon_id)
        overlappingAnnotStart$strand <- myExons2$strand[matchIndex]
        overlappingAnnotStart$exon_start <-
            myExons2$start[matchIndex]
        overlappingAnnotStart$exon_end <- myExons2$end[matchIndex]
        overlappingAnnotStart$exon_cumsum <-
            myExons2$cumSum[matchIndex]

        ### Annoate with transcript coordinats
        overlappingAnnotStartPlus <-
            overlappingAnnotStart[which(overlappingAnnotStart$strand == '+'), ]
        overlappingAnnotStartPlus$transcriptStart <-
            overlappingAnnotStartPlus$exon_cumsum + (
                overlappingAnnotStartPlus$cdsGenomicStart -
                    overlappingAnnotStartPlus$exon_start
            ) + 2

        overlappingAnnotStartMinus <-
            overlappingAnnotStart[which(overlappingAnnotStart$strand == '-'), ]
        overlappingAnnotStartMinus$transcriptStart <-
            overlappingAnnotStartMinus$exon_cumsum + (
                overlappingAnnotStartMinus$exon_end -
                    overlappingAnnotStartMinus$cdsGenomicStart
            ) + 1

        overlappingAnnotStart2 <-
            rbind(overlappingAnnotStartPlus,
                  overlappingAnnotStartMinus)
        #overlappingAnnotStart2 <-
        #    overlappingAnnotStart2[order(
        #        overlappingAnnotStart2$isoform_id,
        #        overlappingAnnotStart2$exon_start,
        #        overlappingAnnotStart2$exon_end
        #    ), ]

        ### Add back in those not overlapping with CDSs      # KVS Dec 2020
        if( orfMethod == 'longest.AnnotatedWhenPossible') {
            allIso <- unique(switchAnalyzeRlist$isoformFeatures$isoform_id)
            overlappingAnnotStart2 <- overlappingAnnotStart2[match(
                allIso, overlappingAnnotStart2$isoform_id
            ),]
            overlappingAnnotStart2$isoform_id <- allIso
        }


        # Update number of steps performed
        nrStepsPerformed <- nrStepsPerformed + 1
    }

    ### Extract nucleotide sequence of the transcripts
    if (!quiet) {
        message(
            paste(
                'Step',
                nrStepsPerformed + 1,
                'of',
                nrStepsToPerform,
                ': Extracting transcript sequences...',
                sep = ' '
            )
        )
    }
    if (TRUE) {
        if( ntAlreadyInSwitchList ) {
            transcriptSequencesDNAstring <- switchAnalyzeRlist$ntSequence

        } else {
            tmpSwitchAnalyzeRlist <- switchAnalyzeRlist

            ### Subset to those analyzed (if annoation is used)
            if (useAnnoated) {
                tmpSwitchAnalyzeRlist$isoform_feature <-
                    tmpSwitchAnalyzeRlist$isoform_feature[which(
                        tmpSwitchAnalyzeRlist$isoform_feature$isoform_id %in%
                            overlappingAnnotStart2$isoform_id
                    ),]
                tmpSwitchAnalyzeRlist$exons <-
                    tmpSwitchAnalyzeRlist$exons[which(
                        tmpSwitchAnalyzeRlist$exons$isoform_id %in%
                            overlappingAnnotStart2$isoform_id
                    ),]
            }

            transcriptSequencesDNAstring <-
                suppressMessages(
                    extractSequence(
                        switchAnalyzeRlist = tmpSwitchAnalyzeRlist,
                        genomeObject = genomeObject,
                        onlySwitchingGenes = FALSE,
                        extractNTseq = TRUE,
                        extractAAseq = FALSE,
                        removeLongAAseq = FALSE,
                        addToSwitchAnalyzeRlist = TRUE,
                        writeToFile = FALSE,
                        quiet = TRUE
                    )$ntSequence
                )

        }

        nrStepsPerformed <- nrStepsPerformed + 1
    }

    ### For each nucleotide sequence identify position of longest ORF with method selected
    if (!quiet) {
        message(
            paste(
                'Step',
                nrStepsPerformed + 1,
                'of',
                nrStepsToPerform,
                ': Locating potential ORFs...',
                sep = ' '
            )
        )
    }
    if (TRUE) {
        if (useAnnoated) {
            overlappingAnnotStartList <-
                split(
                    overlappingAnnotStart2[, c('isoform_id', 'transcriptStart')],
                    f = overlappingAnnotStart2$isoform_id
                )
        } else {
            overlappingAnnotStartList <-
                split(
                    names(transcriptSequencesDNAstring),
                    names(transcriptSequencesDNAstring)
                )
        }

        # Make logics
        useLongest <- orfMethod %in% c('longestAnnotated','longest','longest.AnnotatedWhenPossible')

        ### Find the disired ORF
        transcriptORFs <-
            plyr::llply(
                overlappingAnnotStartList,
                .progress = progressBar,
                function(
                    annoationInfo
                ) { # annoationInfo <- overlappingAnnotStartList[[2]]

                # Extract wanted ORF
                if (useAnnoated) {
                    correspondingSequence <-
                        transcriptSequencesDNAstring[[
                            annoationInfo$isoform_id[1]
                            ]]

                    localORFs <-
                        myAllFindORFsinSeq(
                            dnaSequence = correspondingSequence,
                            codonAnnotation = myCodonDf,
                            filterForPostitions = annoationInfo$transcriptStart
                        )

                    # subset by length
                    localORFs <-
                        localORFs[which(localORFs$length >= minORFlength), ]

                    if (useLongest) {
                        # Longest annoated ORF
                        myMaxORF <-
                            localORFs[which.max(localORFs$length), ]
                    } else {
                        # most upresteam annoated ORF
                        myMaxORF <-
                            localORFs[which.min(localORFs$start), ]
                    }

                } else {
                    correspondingSequence <-
                        transcriptSequencesDNAstring[[annoationInfo]]

                    # longest ORF
                    localORFs <-
                        myAllFindORFsinSeq(
                            dnaSequence = correspondingSequence,
                            codonAnnotation = myCodonDf
                        )

                    # subset by length
                    localORFs <-
                        localORFs[which(localORFs$length >= minORFlength), ]

                    if (useLongest) {
                        # longest ORF
                        myMaxORF <-
                            localORFs[which.max(localORFs$length), ]
                    } else {
                        # most upstream ORF
                        myMaxORF <-
                            localORFs[which.min(localORFs$start), ]
                    }
                }

                # Sanity check
                if (nrow(myMaxORF) == 0) {
                    return(data.frame(
                        start = 1,
                        end = NA,
                        length = 0
                    ))
                } # by having length 0 it will be removed later

                return(myMaxORF)
            })

        myTranscriptORFdf <-
            myListToDf(transcriptORFs, addOrignAsColumn = TRUE)

        nrStepsPerformed <- nrStepsPerformed + 1
    }

    ### Use the obtained ORF coordinats to predict PTC
    if (!quiet) {
        message(
            paste(
                'Step',
                nrStepsPerformed + 1,
                'of',
                nrStepsToPerform,
                ': Scanning for PTCs...',
                sep = ' '
            )
        )
    }
    if (TRUE) {
        ### Extract exon structure for each transcript
        myExons <- as.data.frame(switchAnalyzeRlist$exons[which(
            switchAnalyzeRlist$exons$isoform_id %in%
                switchAnalyzeRlist$isoformFeatures$isoform_id
        ),])
        myExons <- myExons[which(myExons$strand %in% c('+', '-')), ]
        myExons <-
            myExons[which(myExons$isoform_id %in% names(transcriptORFs)), ]
        myExonsSplit <- split(myExons, f = myExons$isoform_id)

        # Loop over all isoforms and extract info
        allIsoforms <-
            split(names(myExonsSplit), names(myExonsSplit))
        ptcResult <-
            plyr::llply(
                allIsoforms,
                .fun = function(isoformName) {
                    # isoformName <- allIsoforms[[1]]
                    # Extract ORF info
                    orfInfo <- transcriptORFs[[isoformName]]
                    if (orfInfo$length == 0) {
                        return(NULL)
                    }

                    # Extract exon info
                    exonInfo <- myExonsSplit[[isoformName]]

                    # do PTC analysis
                    localPTCresult <-
                        analyzeORFforPTC(
                            aORF = orfInfo,
                            exonStructure = exonInfo,
                            PTCDistance = PTCDistance
                        )

                    return(localPTCresult)
                }
            )
        # remove empty once
        ptcResult <-
            ptcResult[which(sapply(ptcResult, function(x)
                ! is.null(x)))]
        if (length(ptcResult) == 0) {
            warning('No ORFs (passing the filtering) were found.\nReturning unmodified switchAnalyzeRlist')
            return(switchAnalyzeRlist)
        }

        myPTCresults <-
            myListToDf(ptcResult, addOrignAsColumn = TRUE)
    }

    ### Add result to switchAnalyzeRlist
    if (TRUE) {
        # merge ORF and PTC analysis together
        myResultDf <-
            dplyr::inner_join(
                myTranscriptORFdf[,c('orign','start','end','length')],
                myPTCresults,
                by = 'orign'
            )
        colnames(myResultDf)[1] <- 'isoform_id'
        colnames(myResultDf)[2:4] <-
            paste(
                'orfTranscipt',
                startCapitalLetter(colnames(myResultDf)[2:4]),
                sep =''
            )

        ### Add NAs
        myResultDf <- merge(
            unique(switchAnalyzeRlist$isoformFeatures[, 'isoform_id', drop = FALSE]),
            myResultDf,
            all.x = TRUE
        )

        # Annotate ORF origin
        myResultDf$orf_origin <- 'Predicted'

        # myResultDf$orfStarExon <- NULL
        # myResultDf$orfEndExon <- NULL

        ### Add result to switch list
        switchAnalyzeRlist$orfAnalysis <- myResultDf

        switchAnalyzeRlist$isoformFeatures$PTC		         <-
            myResultDf$PTC [match(switchAnalyzeRlist$isoformFeatures$isoform_id,
                                  myResultDf$isoform_id)]

        ### Add NT sequences
        switchAnalyzeRlist$ntSequence <- transcriptSequencesDNAstring[which(
            names(transcriptSequencesDNAstring) %in%
                switchAnalyzeRlist$isoformFeatures$isoform_id
        )]

        if (!quiet) {
            message(
                sum(myResultDf$orfStartGenomic != -1, na.rm = TRUE) ,
                " putative ORFs were identified, analyzed and added.",
                sep = ""
            )
        }
    }

    if (!quiet) {
        message('Done')
    }
    return(switchAnalyzeRlist)
}

extractSequence <- function(
    switchAnalyzeRlist,
    genomeObject = NULL,
    onlySwitchingGenes = TRUE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    extractNTseq = TRUE,
    extractAAseq = TRUE,
    removeShortAAseq = TRUE,
    removeLongAAseq  = FALSE,
    alsoSplitFastaFile = FALSE,
    removeORFwithStop = TRUE,
    addToSwitchAnalyzeRlist = TRUE,
    writeToFile = TRUE,
    pathToOutput = getwd(),
    outputPrefix = 'isoformSwitchAnalyzeR_isoform',
    forceReExtraction = FALSE,
    quiet = FALSE
) {
    ### Determine sequence status
    if(TRUE) {
        ntAlreadyInSwitchList <- ! is.null(switchAnalyzeRlist$ntSequence)
        aaAlreadyInSwitchList <- ! is.null(switchAnalyzeRlist$aaSequence)

        if(forceReExtraction) {
            ntAlreadyInSwitchList <- FALSE
            aaAlreadyInSwitchList <- FALSE
        }
    }

    ### Check input
    if (TRUE) {
        # Test input data class
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(
                'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
            )
        }
        if( ! ntAlreadyInSwitchList ) {
            if (class(genomeObject) != 'BSgenome') {
                stop('The genomeObject argument must be a BSgenome')
            }
        }

        # Test what to extract
        if (!extractNTseq & !extractAAseq) {
            stop(
                'At least one of \'extractNTseq\' or \'extractNTseq\' must be true (else using this function have no purpose)'
            )
        }

        # How to repport result
        if (!addToSwitchAnalyzeRlist & !writeToFile) {
            stop(
                'At least one of \'addToSwitchAnalyzeRlist\' or \'writeToFile\' must be true (else this function outputs nothing)'
            )
        }

        # Are ORF annotated
        if (extractAAseq) {
            if (!'orfAnalysis' %in% names(switchAnalyzeRlist)) {
                stop('Please run the \'addORFfromGTF()\' (and if nessesary \'analyzeNovelIsoformORF()\') function(s) to detect ORFs')
            }

            if( ! is.null(switchAnalyzeRlist$orfAnalysis$orf_origin) ) {
                if ( any( switchAnalyzeRlist$orfAnalysis$orf_origin == 'not_annotated_yet' )) {
                    stop('Some ORFs have not been annotated yet. Please run analyzeNovelIsoformORF() and try again.')
                }
            }

        }

        # If switches are annotated
        if (onlySwitchingGenes) {
            if (all(is.na(
                switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
            )) &
            all(is.na(
                switchAnalyzeRlist$isoformFeatures$gene_switch_q_value
            ))) {
                stop(
                    'If only switching genes should be outputted please run the \'isoformSwitchTestDEXSeq\' or \'isoformSwitchTestDRIMSeq\' function first and try again'
                )
            }
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

        if( !is.logical(alsoSplitFastaFile) ){
            stop('The \'alsoSplitFastaFile\' argument must be either TRUE or FALSE')
        }
        if( alsoSplitFastaFile ) {
            if( ! removeShortAAseq ) {
                warning('Since you are using the alsoSplitFastaFile you probably also want to use the \'removeShortAAseq\' option.')
            }
            if( ! removeLongAAseq ) {
                warning('Since you are using the alsoSplitFastaFile you probably also want to use the \'removeLongAAseq\' option.')
            }
        }

        if( !is.logical(forceReExtraction)) {
            stop('\'forceReExtraction\' must be either TRUE or FALSE')
        }

        nrAnalysisToMake <- 2 + as.integer(extractAAseq)
        startOfAnalysis <- 1
    }

    ### Extract NT sequence (needed for AA extraction so always nessesary)
    if (TRUE) {
        if (!quiet) {
            message(
                paste(
                    "Step",
                    startOfAnalysis ,
                    "of",
                    nrAnalysisToMake,
                    ": Extracting transcript nucleotide sequences...",
                    sep = " "
                )
            )
        }

        # extract switching isoforms
        if (onlySwitchingGenes) {
            isoResTest <-
                any(!is.na(
                    switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
                ))
            if (isoResTest) {
                switchingGenes <-
                    unique(switchAnalyzeRlist$isoformFeatures$gene_id [which(
                        switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value <
                            alpha &
                            abs(switchAnalyzeRlist$isoformFeatures$dIF) >
                            dIFcutoff
                    )])
            } else {
                switchingGenes <-
                    unique(switchAnalyzeRlist$isoformFeatures$gene_id [which(
                        switchAnalyzeRlist$isoformFeatures$gene_switch_q_value <
                            alpha     &
                            abs(switchAnalyzeRlist$isoformFeatures$dIF) >
                            dIFcutoff
                    )])
            }
            if (length(switchingGenes) == 0) {
                stop(
                    'No switching genes were found. Pleasae turn off \'onlySwitchingGenes\' and try again.'
                )
            }
            switchingIsoforms <-
                unique(switchAnalyzeRlist$isoformFeatures$isoform_id[which(
                    switchAnalyzeRlist$isoformFeatures$gene_id %in%
                        switchingGenes
                )])
        }

        # Extract nuclotide sequences
        if(   ntAlreadyInSwitchList ) {
            if (onlySwitchingGenes) {
                transcriptSequencesDNAstring <- switchAnalyzeRlist$ntSequence[which(
                    names(switchAnalyzeRlist$ntSequence) %in% switchingIsoforms
                )]
            } else {
                transcriptSequencesDNAstring <- switchAnalyzeRlist$ntSequence
            }

        }
        if( ! ntAlreadyInSwitchList ) {
            ### Extract exon GRanges
            if (onlySwitchingGenes) {
                # only extract transcript info for the switching genes
                myExonGranges <-
                    switchAnalyzeRlist$exons[which(
                        switchAnalyzeRlist$exons$isoform_id %in% switchingIsoforms
                    ),]
            } else {
                myExonGranges <- switchAnalyzeRlist$exons[which(
                    switchAnalyzeRlist$exons$isoform_id %in%
                        switchAnalyzeRlist$isoformFeatures$isoform_id
                ),]
            }

            myExonGranges <-
                myExonGranges[which(
                    as.character(myExonGranges@strand) %in% c('+', '-')
                ) , ]

            # update grange levels (migth be nessesary if it is a subset)
            seqlevels(myExonGranges) <-
                unique( as.character(seqnames(myExonGranges)@values) ) # nessesary to make sure seqlevels not pressented in the input data casuses problems - if presented with a subset for example


            ### Check whether isoform annotation and genome fits together
            genomeSeqs <- seqnames(genomeObject)
            grSeqs <- seqlevels(myExonGranges)

            # check overlap
            chrOverlap <- intersect(grSeqs, genomeSeqs)

            if (length(chrOverlap) == 0) {
                ### Try correct chr names of GRanges
                if (any(!grepl('^chr', seqlevels(myExonGranges)))) {
                    # Correct all chromosomes
                    seqlevels(myExonGranges) <-
                        unique(paste("chr", seqnames(myExonGranges), sep = ""))
                    # Correct Mitochondria
                    seqlevels(myExonGranges) <-
                        sub('chrMT', 'chrM', seqlevels(myExonGranges))

                    # check overlap again
                    grSeqs <- seqlevels(myExonGranges)
                    chrOverlap <- intersect(grSeqs, genomeSeqs)

                    # Correct chr names genome
                } else if (any(!grepl('^chr', genomeSeqs))) {
                    # Correct all chromosomes
                    seqnames(genomeObject) <-
                        paste("chr", seqnames(genomeObject), sep = "")
                    # Correct Mitochondria
                    seqnames(genomeObject) <-
                        sub('chrMT', 'chrM', seqnames(genomeObject))

                    # check overlap again
                    genomeSeqs <- seqnames(genomeObject)
                    chrOverlap <- intersect(grSeqs, genomeSeqs)

                }

                if (length(chrOverlap) == 0) {
                    stop(
                        'The genome supplied to genomeObject have no seqnames in common with the genes in the switchAnalyzeRlist'
                    )
                }
            }

            ### remove those transcripts that does not have a corresponding chromosme
            notOverlappingIndex <- which(!grSeqs %in% genomeSeqs)
            if (length(notOverlappingIndex)) {
                toRemove <-
                    which(seqnames(myExonGranges) %in% grSeqs[notOverlappingIndex])

                if (length(toRemove)) {
                    nrTranscripts <-
                        length(unique(myExonGranges$isoform_id[toRemove]))
                    warning(
                        paste(
                            nrTranscripts,
                            'transcripts was removed due to them being annotated to chromosomes not found in the refrence genome object',
                            sep = ' '
                        )
                    )

                    myExonGranges <- myExonGranges[-toRemove , ]
                }
            }

            # extract exon sequence and split into transcripts
            myExonSequences <- getSeq(genomeObject, myExonGranges)
            names(myExonSequences) <- myExonGranges$isoform_id

            ### In a strand specific manner combine exon sequenes (strand specific is nessesry because they should be reversed - else they are combined in the '+' order)
            myStrandData <-
                unique(
                    data.frame(
                        isoform_id = myExonGranges$isoform_id,
                        strand = as.character(myExonGranges@strand),
                        stringsAsFactors = FALSE
                    )
                )
            plusStrandTransripts <-
                myStrandData$isoform_id[which(myStrandData$strand == '+')]
            minusStrandTransripts <-
                myStrandData$isoform_id[which(myStrandData$strand == '-')]

            # Collaps exons of plus strand transcripts
            myPlusExonSequences <-
                myExonSequences[which(
                    names(myExonSequences) %in% plusStrandTransripts
                ), ]
            myPlusExonSequences <-
                split(myPlusExonSequences, f = names(myPlusExonSequences))
            myPlusExonSequences <-
                lapply(myPlusExonSequences, unlist) # does not work in the R package

            # Collaps exons of minus strand transcripts
            myMinusExonSequences <-
                myExonSequences[which(
                    names(myExonSequences) %in% minusStrandTransripts
                ), ]
            myMinusExonSequences <- rev(myMinusExonSequences)
            myMinusExonSequences <-
                split(myMinusExonSequences, f = names(myMinusExonSequences))
            myMinusExonSequences <- lapply(myMinusExonSequences, unlist)

            # combine strands
            transcriptSequencesDNAstring <-
                DNAStringSet(c(myPlusExonSequences, myMinusExonSequences))

        }

        startOfAnalysis <- startOfAnalysis + 1
    }

    ### Extract protein sequence of the identified ORFs
    if (extractAAseq) {
        if (!quiet) {
            message(
                paste(
                    "Step",
                    startOfAnalysis ,
                    "of",
                    nrAnalysisToMake,
                    ": Extracting ORF AA sequences...",
                    sep = " "
                )
            )
        }

        ### Extract switchAnalyzeRlist ORF annotation and filter for those I need
        if(TRUE) {
            switchORFannotation <-
                unique(data.frame(switchAnalyzeRlist$orfAnalysis[, c(
                    'isoform_id',
                    "orfTransciptStart",
                    'orfTransciptEnd',
                    "orfTransciptLength",
                    "PTC"
                )]))
            switchORFannotation <-
                switchORFannotation[which(!is.na(switchORFannotation$PTC)), ]
            switchORFannotation <-
                switchORFannotation[which(
                    switchORFannotation$isoform_id %in%
                        names(transcriptSequencesDNAstring)), ]
            switchORFannotation <-
                switchORFannotation[which(
                    switchORFannotation$orfTransciptStart != 0
                ), ]

        }

        ### Extract AA sequences
        if(   aaAlreadyInSwitchList ) {
            transcriptORFaaSeq <- switchAnalyzeRlist$aaSequence[which(
                names(switchAnalyzeRlist$aaSequence) %in% names(transcriptSequencesDNAstring)
            )]
        }
        if( ! aaAlreadyInSwitchList ) {
            ### Reorder transcript sequences
            transcriptSequencesDNAstringInData <-
                transcriptSequencesDNAstring[na.omit(match(
                    x = switchORFannotation$isoform_id,
                    table = names(transcriptSequencesDNAstring)
                )), ]
            if (!all(
                names(transcriptSequencesDNAstringInData) ==
                switchORFannotation$isoform_id
            )) {
                stop('Somthing went wrong in sequence extraction - contract developer')
            }

            ### Test whether the annotation agrees
            switchORFannotation$lengthOK <-
                switchORFannotation$orfTransciptEnd <=
                width(transcriptSequencesDNAstringInData)[match(
                    switchORFannotation$isoform_id,
                    names(transcriptSequencesDNAstringInData)
                )]
            if (any(!switchORFannotation$lengthOK)) {
                warning(
                    paste(
                        'There were',
                        sum(!switchORFannotation$lengthOK),
                        'cases where the annotated ORF were longer than the exons annoated - these cases will be ommitted'
                    )
                )

                # Subset data
                switchORFannotation <-
                    switchORFannotation[which(switchORFannotation$lengthOK), ]
                transcriptSequencesDNAstringInData <-
                    transcriptSequencesDNAstringInData[na.omit(match(
                        x = switchORFannotation$isoform_id,
                        table = names(transcriptSequencesDNAstringInData)
                    )), ]
            }

            ### Get corresponding protein sequence
            # Use the predicted ORF coordinats to extract the nt sequence of the ORF
            transcriptORFntSeq <-
                XVector::subseq(
                    transcriptSequencesDNAstringInData,
                    start = switchORFannotation$orfTransciptStart,
                    width = switchORFannotation$orfTransciptLength
                )

            # translate ORF nucleotide to aa sequence
            transcriptORFaaSeq <-
                suppressWarnings(
                    Biostrings::translate(x = transcriptORFntSeq, if.fuzzy.codon = 'solve')
                ) # supress warning is nessesary because isoformSwitchAnalyzeR allows ORFs to exceed the transcript - which are by default ignored and just gives a warning

            ### Trim (potential) last stop codon
            transcriptORFaaSeq <- Biostrings::trimLRPatterns(
                Lpattern = '',
                Rpattern = "*",
                subject = transcriptORFaaSeq
            )

            ### Check ORFs for stop codons
            stopData <- data.frame(
                isoform_id = names(transcriptORFaaSeq),
                stopCodon = Biostrings::vcountPattern(pattern = '*', transcriptORFaaSeq),
                stringsAsFactors = FALSE
            )
            stopDataToRemove <- stopData[which(stopData$stopCodon > 0), ]

            if (nrow(stopDataToRemove)) {
                if (removeORFwithStop) {
                    warning(
                        paste(
                            'There were',
                            nrow(stopDataToRemove),
                            'isoforms where the amino acid sequence had a stop codon before the annotated stop codon. These was removed.',
                            sep = ' '
                        )
                    )

                    ### Remove PTC annotation
                    switchAnalyzeRlist$isoformFeatures$PTC[which(
                        switchAnalyzeRlist$isoformFeatures$isoform_id %in% stopDataToRemove$isoform_id
                    )] <- NA

                    ### Remove ORF annoation
                    switchAnalyzeRlist$orfAnalysis[which(
                        switchAnalyzeRlist$orfAnalysis$isoform_id %in% stopDataToRemove$isoform_id
                    ), which( ! colnames(switchAnalyzeRlist$orfAnalysis) %in% c('isoform_id','orf_origin')) ] <- NA

                    ### Remove sequence
                    transcriptORFaaSeq <-
                        transcriptORFaaSeq[which(!names(transcriptORFaaSeq) %in% stopDataToRemove$isoform_id)]

                    switchORFannotation <- switchORFannotation[which(
                        ! switchORFannotation$isoform_id %in% stopDataToRemove$isoform_id
                    ),]

                } else {
                    warning(
                        paste(
                            'There were',
                            nrow(stopDataToRemove),
                            'isoforms where the amino acid sequence had a stop codon before the annotated stop codon. These was NOT removed in accodance with the \'removeORFwithStop\' argument.',
                            sep = ' '
                        )
                    )
                }
            }



        }

        startOfAnalysis <- startOfAnalysis + 1

    }


    ### If enabled write fasta file(s) (after filteri g)
    seqWasTrimmed <- FALSE
    if (writeToFile) {
        if (!quiet) {
            message(
                paste(
                    "Step",
                    startOfAnalysis ,
                    "of",
                    nrAnalysisToMake,
                    ": Preparing output...",
                    sep = " "
                )
            )
        }

        ### add / if directory
        if (file.exists(pathToOutput)) {
            pathToOutput <- paste(pathToOutput, '/', sep = '')
        } else {
            stop('The path supplied to \'pathToOutput\' does not seem to exist')
        }

        # Nucleotides
        if (extractNTseq) {
            writeXStringSet(
                transcriptSequencesDNAstring,
                filepath = paste(
                    pathToOutput,
                    outputPrefix,
                    '_nt.fasta',
                    sep = ''
                ),
                format = 'fasta'
            )
        }

        # Amino Acids
        if (extractAAseq) {

            ### Filter if nessesary
            if (removeLongAAseq | removeShortAAseq) {

                switchORFannotation$toBeTrimmed <- FALSE

                ### remove to long sequences
                if(   removeLongAAseq ) {
                    ### Filter lengths
                    switchORFannotation$aaLength <- width(transcriptORFaaSeq)
                    switchORFannotation$toBeTrimmed <- switchORFannotation$aaLength > 1000


                    ### Make new lengths
                    switchORFannotation$newAAlength <- switchORFannotation$aaLength
                    switchORFannotation$newAAlength[which(
                        switchORFannotation$toBeTrimmed
                    )] <- 1000   # <- EBI's current limmit

                    ### Annotate the trimming
                    if(any( switchORFannotation$toBeTrimmed )) {
                        seqWasTrimmed <- TRUE


                        trimmedCases <- switchORFannotation[which(
                            switchORFannotation$toBeTrimmed
                        ),]

                        ### Calculate genomic positions of trimmed regions
                        if(TRUE) {
                            trimmedCases$trimmedTranscriptStart <-
                                (1001  * 3 - 2) + trimmedCases$orfTransciptStart - 1
                            trimmedCases$trimmedTranscriptEnd <- trimmedCases$orfTransciptEnd

                            ### convert from transcript to genomic coordinats
                            # extract exon data
                            myExons <-
                                as.data.frame(switchAnalyzeRlist$exons[which(
                                    switchAnalyzeRlist$exons$isoform_id %in% trimmedCases$isoform_id
                                ), ])
                            myExonsSplit <- split(myExons, f = myExons$isoform_id)

                            # loop over the individual transcripts and extract the genomic coordiants of the domain and also for the active residues (takes 2 min for 17000 rows)
                            trimmedCases <-
                                plyr::ddply(
                                    trimmedCases[,c('isoform_id','newAAlength','trimmedTranscriptStart','trimmedTranscriptEnd')],
                                    .progress = 'none',
                                    .variables = 'isoform_id',
                                    .fun = function(aDF) {
                                        # aDF <- trimmedCases[1,]

                                        transcriptId <- aDF$isoform_id[1]
                                        localExons <-
                                            as.data.frame(myExonsSplit[[transcriptId]])

                                        # extract domain allignement
                                        localORFalignment <- aDF
                                        colnames(localORFalignment)[match(
                                            x = c('trimmedTranscriptStart', 'trimmedTranscriptEnd'),
                                            table = colnames(localORFalignment)
                                        )] <- c('start', 'end')

                                        # loop over domain alignment (migh be several)
                                        orfPosList <- list()
                                        for (j in 1:nrow(localORFalignment)) {
                                            domainInfo <-
                                                convertCoordinatsTranscriptToGenomic(
                                                    transcriptCoordinats =  localORFalignment[j, ],
                                                    exonStructure = localExons
                                                )

                                            orfPosList[[as.character(j)]] <- domainInfo

                                        }
                                        orfPosDf <- do.call(rbind, orfPosList)

                                        return(cbind(aDF, orfPosDf))
                                    }
                                )

                        }

                    }

                    ### Trim
                    transcriptORFaaSeq2 <- XVector::subseq(
                        transcriptORFaaSeq,
                        start = 1,
                        width = switchORFannotation$newAAlength
                    )
                }
                if( ! removeLongAAseq) {
                    transcriptORFaaSeq2 <- transcriptORFaaSeq
                }

                ### Remove to short sequences
                if( removeShortAAseq ) {
                    transcriptORFaaSeq2 <-
                        transcriptORFaaSeq2[which(
                            width(transcriptORFaaSeq2) > 10
                        )]
                }

                message(
                    paste(
                        'The \'removeLongAAseq\' and \'removeShortAAseq\' arguments:\n',
                        'Removed :',
                        length(transcriptORFaaSeq) - length(transcriptORFaaSeq2),
                        'isoforms.\n',
                        'Trimmed :',
                        sum(switchORFannotation$toBeTrimmed),
                        'isoforms (to only contain the first 1000 AA)',
                        sep=' '
                    )
                )
            } else {
                transcriptORFaaSeq2 <- transcriptORFaaSeq
            }

            ### Write file(s)
            if(   alsoSplitFastaFile ) {
                ### Make index
                l <- length(transcriptORFaaSeq2)

                maxfileSizes <- 500
                nFiles <- ceiling(l / maxfileSizes)
                seqWithinEachFile <-  ceiling(l / nFiles)

                indexVec <- unique( c( seq(
                    from = 1,
                    to = l,
                    by = seqWithinEachFile # Max in PFAM Jan 2019
                ), l))

                indexDf <- data.frame(
                    start = indexVec[-length(indexVec)],
                    end = indexVec[-1]
                )
                n <- nrow(indexDf)
                indexDf$file <- paste0('_subset_', 1:n,'_of_',n)

                ### loop over index and make files
                tmp <- plyr::ddply(indexDf, .variables = 'file', .fun = function(aDF) { # aDF <- indexDf[1,]
                    writeXStringSet(
                        transcriptORFaaSeq2[ aDF$start:aDF$end ],
                        filepath = paste(
                            pathToOutput,
                            outputPrefix,
                            '_AA',
                            aDF$file,
                            '.fasta',
                            sep = ''
                        ),
                        format = 'fasta'
                    )
                })
                message(paste(
                    'The \'alsoSplitFastaFile\' caused',
                    nrow(indexDf),
                    'fasta files, each with a subset of the data, to be created (each named X of Y).'
                ))

                ### ALso write full
                writeXStringSet(
                    transcriptORFaaSeq2,
                    filepath = paste(
                        pathToOutput,
                        outputPrefix,
                        '_AA_complete.fasta',
                        sep = ''
                    ),
                    format = 'fasta'
                )
            }
            if( ! alsoSplitFastaFile ) {
                ### Write full
                writeXStringSet(
                    transcriptORFaaSeq2,
                    filepath = paste(
                        pathToOutput,
                        outputPrefix,
                        '_AA.fasta',
                        sep = ''
                    ),
                    format = 'fasta'
                )
            }

        }
    }

    if (!quiet) {
        message('Done')
    }

    ### Add sequences to switchAnalyzeRlist
    if (addToSwitchAnalyzeRlist) {
        if (extractNTseq) {
            switchAnalyzeRlist$ntSequence <- transcriptSequencesDNAstring
        }
        if (extractAAseq) {
            switchAnalyzeRlist$aaSequence <- transcriptORFaaSeq
        }

    }

    ### Add run info
    if( TRUE ) {
        if( is.null(switchAnalyzeRlist$runInfo) ) {
            switchAnalyzeRlist$runInfo <- list()
        }

        if(seqWasTrimmed) {
            switchAnalyzeRlist$orfAnalysis$wasTrimmed <- switchORFannotation$toBeTrimmed[match(
                switchAnalyzeRlist$orfAnalysis$isoform_id, switchORFannotation$isoform_id
            )]

            switchAnalyzeRlist$orfAnalysis$trimmedStartGenomic <- trimmedCases$pfamStartGenomic[match(
                switchAnalyzeRlist$orfAnalysis$isoform_id, trimmedCases$isoform_id
            )]
        }

        switchAnalyzeRlist$runInfo$extractSequence <- list(
            removeShortAAseq = removeShortAAseq,
            removeLongAAseq = seqWasTrimmed
        )

    }
    return(switchAnalyzeRlist)
}


addORFfromGTF <- function(
    ### Core arguments
    switchAnalyzeRlist,
    pathToGTF,

    ### Advanced argument
    overwriteExistingORF = FALSE,
    onlyConsiderFullORF = FALSE,
    removeNonConvensionalChr = FALSE,
    ignoreAfterBar = TRUE,
    ignoreAfterSpace = TRUE,
    ignoreAfterPeriod = FALSE,
    PTCDistance = 50,
    quiet = FALSE
) {
    ### Improvements:
    # 1) have importRdata() save ignoreAfterBar, ignoreAfterSpace and ignoreAfterPeriod options and re-use
    # 2) Test overlap. But not untill 1) is done

    ### Test input
    if(TRUE) {
        # Test input data class
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(
                'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
            )
        }

        if( !is.null(switchAnalyzeRlist$orfAnalysis ) ) {
            if( ! all(is.na(switchAnalyzeRlist$orfAnalysis$orfTransciptStart ))) {
                if( ! overwriteExistingORF) {
                    stop('ORF appears to already be annotated. Set \'overwriteExistingORF=TRUE\' to overwrite anyway.')
                }
            }
        }

    }

    ### Import GTF and CDS
    if(TRUE) {
        if (!quiet) {
            message('Step 1 of 2: importing GTF (this may take a while)...')
        }

        suppressWarnings(
            gtfSwichList <- importGTF(
                pathToGTF = pathToGTF,
                addAnnotatedORFs = TRUE,
                onlyConsiderFullORF = onlyConsiderFullORF,
                removeNonConvensionalChr = FALSE,
                ignoreAfterBar = ignoreAfterBar,
                ignoreAfterSpace = ignoreAfterSpace,
                ignoreAfterPeriod = ignoreAfterPeriod,
                removeTECgenes = FALSE,
                PTCDistance = PTCDistance,
                removeFusionTranscripts = FALSE,
                quiet = TRUE
            )
        )

        ### Test ORF annotaion
        if( is.null(gtfSwichList$orfAnalysis) ) {
            stop('No CDS information was extracted from the GTF file.')
        }
    }

    ### Add CDS to switchList
    if(TRUE) {
        if (!quiet) {
            message('Step 2 of 2: Adding ORF...')
        }

        ### Add to switch list
        allIso <- unique(switchAnalyzeRlist$isoformFeatures$isoform_id)

        switchAnalyzeRlist$orfAnalysis <- gtfSwichList$orfAnalysis[match(
            allIso, gtfSwichList$orfAnalysis$isoform_id
        ),]
        switchAnalyzeRlist$orfAnalysis$isoform_id <- allIso

        switchAnalyzeRlist$isoformFeatures$PTC <-
            switchAnalyzeRlist$orfAnalysis$PTC[match(
                switchAnalyzeRlist$isoformFeatures$isoform_id,
                switchAnalyzeRlist$orfAnalysis$isoform_id
            )]

        ### Add origin
        switchAnalyzeRlist$orfAnalysis$orf_origin[which(
            is.na(switchAnalyzeRlist$orfAnalysis$orf_origin)
        )] <- 'not_annotated_yet'

        ### Extract summary numbers
        if(TRUE) {
            ### Extract info
            inSl <- unique(switchAnalyzeRlist$isoformFeatures$isoform_id)
            inORF <- unique(switchAnalyzeRlist$orfAnalysis$isoform_id[which(
                switchAnalyzeRlist$orfAnalysis$orf_origin != 'not_annotated_yet'
            )])
            hasORF <- unique(switchAnalyzeRlist$orfAnalysis$isoform_id[which(
                switchAnalyzeRlist$orfAnalysis$orf_origin != 'not_annotated_yet' &
                ! is.na(switchAnalyzeRlist$orfAnalysis$orfTransciptStart)
            )])
            notAnalysed <- unique(switchAnalyzeRlist$orfAnalysis$isoform_id[which(
                switchAnalyzeRlist$orfAnalysis$orf_origin == 'not_annotated_yet'
            )])

            nInSL <- length(inSl)
            nAdded <- length(inORF)
            nWithOrf <- length(hasORF)

            knownIso <- unique(
                switchAnalyzeRlist$isoformFeatures$isoform_id[which(
                    ! is.na(switchAnalyzeRlist$isoformFeatures$gene_name)
                )]
            )
            knownIsoAdded <- intersect(
                switchAnalyzeRlist$isoformFeatures$isoform_id[which(
                    ! is.na(switchAnalyzeRlist$isoformFeatures$gene_name)
                )],
                inORF
            )
            knownNotAdded <- intersect(knownIso, notAnalysed)

            nNotAnnoateted <- nInSL - nAdded
        }

        ### Repport numbers
        if (!quiet) {
            ### Write message
            message(paste0(
                '    Added ORF info (incl info about isoforms annotated as not having an ORF) to ',
                nAdded, ' isoforms.',
                '\n        This correspond to ', round(nAdded /nInSL * 100, digits=2), '% of isoforms in the switchAnalyzeRlist.',
                '\n            Which includes ', round(length(knownIsoAdded) / length(knownIso) * 100, digits=2), '% of isoforms from annotated genes (novel isoforms not counted) in the switchAnalyzeRlist.'
            ))

            ### Continue message in cases some where not annotated
            if( nNotAnnoateted > 0) {

                isoNotAnnoated <- setdiff(inSl, inORF)

                message(paste0(
                    '    Isoforms with no ORF annoation are either due to incompatible annotation versions or novel isoforms.',
                    '\n        We estimate ', round( (( (length(knownIso) - length(knownNotAdded)) / length(knownIso) )) * 100, digits=2), '% of isoforms not annotated with ORF are from novel genes.',
                    '\n        Examples of isoforms where ORF annoation is still missing are:',
                    '\n            ', paste(sample(isoNotAnnoated, size = min(c(3, nNotAnnoateted))), collapse = ', ')
                ))
            }

        }

        ### Test overlap
        if(TRUE) {
            if( nWithOrf /nInSL == 0) {
                stop(str_c(
                    'No ORFs could be added to the switchAnalyzeRlist.',
                    ' Please ensure GTF file have CDS info and that isoform ids match.'
                ))
            }

            if( nWithOrf /nInSL < 0.5 ) {
                if( nWithOrf /nInSL < 0.1 ) {
                    warning(paste0(
                        'It seems ORF was only added to an unlikely small fraction of isoforms (less than 10%).',
                        'If you are not sure this is on purpose something went wrong and the files are not matching.'
                    ))
                } else {
                    warning(paste0(
                        'It seems ORF was added to an small fraction of isoforms (less than 50%).',
                        'You might want to dobule check the ratio of know/novel transcripts to ensure this is not a problem with files not matching.'
                    ))
                }
            }

        }
    }

    ### Return
    if (!quiet) {
        message('Done.')
        return(switchAnalyzeRlist)
    }
}

analyzeNovelIsoformORF <- function(
    ### Core arguments
    switchAnalyzeRlist,
    analysisAllIsoformsWithoutORF, # also analyse all those annoatated as without CDS in ref annottaion
    genomeObject = NULL,

    ### Advanced argument
    minORFlength = 100,
    orfMethod = 'longest.AnnotatedWhenPossible',
    PTCDistance = 50,
    startCodons = "ATG",
    stopCodons = c("TAA", "TAG", "TGA"),
    showProgress = TRUE,
    quiet = FALSE
) {
    ### Improvements
    #

    ### Test input
    if(TRUE) {
        # Test input data class
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(
                'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
            )
        }

        if( is.null(switchAnalyzeRlist$orfAnalysis) ) {
            stop('No ORF annotation pressent. Run addORFfromGTF() first and try again.')
        }
        nWithout <- sum(switchAnalyzeRlist$orfAnalysis$orf_origin == 'not_annotated_yet')
        nWithout <- nWithout + (sum(!is.na(switchAnalyzeRlist$orfAnalysis$PTC)) * as.integer(analysisAllIsoformsWithoutORF))

        if( nWithout == 0) {
            stop('There appear not to be any isoforms not already annotated with ORFs - meaning there is no need to run this function')
        }

        nWith <- sum(!is.na(switchAnalyzeRlist$orfAnalysis$orfTransciptStart))
        if(nWith == 0) {
            stop('While ORFs seems to have beeen analysed no ORFs were identified. Consider using less strict ORF detection criteria or if something else about ORF detection could have gone wrong.')
        }

        nKnown <- sum(
            switchAnalyzeRlist$orfAnalysis$orf_origin == 'Annotation'
        )
        if(nKnown == 0) {
            stop('There appear to be no known ORFs annotated and hence this function can not be used.')
        }


        ntAlreadyInSwitchList <- ! is.null(switchAnalyzeRlist$ntSequence)
        if( ! ntAlreadyInSwitchList ) {
            if (class(genomeObject) != 'BSgenome') {
                stop('The genomeObject argument must be supplied with a BSgenome object when isoform sequences are not already stored in the switchAnalyzeRlist.')
            }
        }


        # method choice
        if (length(orfMethod) != 1) {
            stop(paste(
                'The \'orfMethod\' must be one of \'mostUpstreamAnnoated\',',
                '\'mostUpstream\', \'longestAnnotated\', \'longest.AnnotatedWhenPossible\' or \'longest\',',
                'not a vector'
            ))
        }
        if (!orfMethod %in% c(
            'mostUpstreamAnnoated',
            'mostUpstream',
            'longest',
            'longestAnnotated',
            'longest.AnnotatedWhenPossible'
        )
        ) {
            stop(paste(
                'The \'orfMethod\' argument must be either of',
                '\'mostUpstreamAnnoated\',\'mostUpstream\',',
                ' \'longestAnnotated\', \'longest.AnnotatedWhenPossible\' or \'longest\' indicating whether',
                'to use the the most upstream annoated start codon or',
                'the longest ORF respectively'
            ))
        }
    }

    ### Step 1 : extact existing ORFs as a CDS object
    if(TRUE) {
        if( orfMethod %in% c(
            'longestAnnotated',
            'mostUpstreamAnnoated',
            'longest.AnnotatedWhenPossible'
        ) ) {
            if (!quiet) {
                message('Step 0 of 4 : Extracting CDS from already annotated isoforms...')
            }

            knownCds <-
                switchAnalyzeRlist$orfAnalysis %>%
                dplyr::filter(
                    orf_origin == 'Annotation',
                    !is.na(orfStartGenomic)
                ) %>%
                dplyr::select(isoform_id, orfStartGenomic, orfEndGenomic)

            knownCds$chrom <- as.character(seqnames(switchAnalyzeRlist$exons[match(
                knownCds$isoform_id, switchAnalyzeRlist$exons$isoform_id
            ),]))
            knownCds$strand <- as.character(strand(switchAnalyzeRlist$exons[match(
                knownCds$isoform_id, switchAnalyzeRlist$exons$isoform_id
            ),]))

            ### Change cds order to genome oriented instead of transcript oriented
            knownCds$cdsStart <- pmin(knownCds$orfStartGenomic, knownCds$orfEndGenomic)
            knownCds$cdsEnd   <- pmax(knownCds$orfStartGenomic, knownCds$orfEndGenomic)

            ### Correct as input is expected to be 0 based
            knownCds$cdsStart <- knownCds$cdsStart - 1

            knownCds <- CDSSet(knownCds)
        } else {
            knownCds <- NULL
        }

    }

    ### Step 2 : use regular analyseORF with the CDS object on novel isoforms
    if(TRUE) {
        ### Extract switchAnalyzeRlist with unannotated isoforms
        nonAnnotatedIso <- switchAnalyzeRlist$orfAnalysis$isoform_id[which(
            switchAnalyzeRlist$orfAnalysis$orf_origin == 'not_annotated_yet'
        )]

        if(analysisAllIsoformsWithoutORF) {
            nonAnnotatedIso <- c(
                nonAnnotatedIso,
                switchAnalyzeRlist$orfAnalysis$isoform_id[which(
                    switchAnalyzeRlist$orfAnalysis$orf_origin == 'Annotation' &
                        is.na(switchAnalyzeRlist$orfAnalysis$orfTransciptStart)
                )]
            )
        }

        unannotatedSl <- subsetSwitchAnalyzeRlist(
            switchAnalyzeRlist = switchAnalyzeRlist,
            switchAnalyzeRlist$isoformFeatures$isoform_id %in% nonAnnotatedIso
        )

        ### Analyse ORFs
        annotatedSl <- analyzeORF(
            switchAnalyzeRlist = unannotatedSl,
            genomeObject = genomeObject,
            orfMethod = orfMethod,
            cds = knownCds,
            minORFlength = minORFlength,
            PTCDistance = PTCDistance,
            startCodons = startCodons,
            stopCodons = stopCodons,
            showProgress = showProgress,
            quiet = quiet
        )

        ### Overwrite analysed results in switchAnalyzeRlist
        switchAnalyzeRlist$orfAnalysis[match(
            annotatedSl$orfAnalysis$isoform_id, switchAnalyzeRlist$orfAnalysis$isoform_id
        ),] <- annotatedSl$orfAnalysis

        switchAnalyzeRlist$isoformFeatures$PTC <-
            switchAnalyzeRlist$orfAnalysis$PTC[match(
                switchAnalyzeRlist$isoformFeatures$isoform_id,
                switchAnalyzeRlist$orfAnalysis$isoform_id
            )]
    }

    return(switchAnalyzeRlist)
}


### Helper functions
myAllFindORFsinSeq <- function(
    dnaSequence, # A DNAString object containing the DNA nucleotide sequence of to analyze for open reading frames
    codonAnnotation = data.frame(
        codons = c("ATG", "TAA", "TAG", "TGA"),
        meaning = c('start', 'stop', 'stop', 'stop'),
        stringsAsFactors = FALSE
    ), # a data.frame with two collums: 1) A collumn called \'codons\' containing a vector of capitalized three-letter strings with codons to analyze. 2) A collumn called \'meaning\' containing the corresponding meaning of the codons. These must be either \'start\' or \'stop\'. See defult data.frame for example. Default are canonical \'ATG\' as start start codons and \'TAA\', \'TAG\', \'TGA\' as stop codons.
    filterForPostitions = NULL # A vector of which transcipt start site potitions to extract

) {
    ### To do
    # might be made more efficeint using matchPDict()
    # Find all codons of interes in the supplied dnaSequence
    myFindPotentialStartsAndStops <- function(dnaSequence, codons) {
        # use matchPattern() to find the codons
        myCodonPositions <-
            sapply(codons, function(x)
                matchPattern(x, dnaSequence))

        # extract the position of the codons
        myCodonPositions <-
            lapply(myCodonPositions, function(x)
                data.frame(position = start(x@ranges)))
        myCodonPositions <-
            myListToDf(
                myCodonPositions,
                ignoreColNames = TRUE,
                addOrignAsRowNames = FALSE,
                addOrignAsColumn = TRUE,
                addOrgRownames = FALSE
            )

        # Massage the data.frame
        colnames(myCodonPositions)[2] <- 'codon'
        myCodonPositions <-
            myCodonPositions[sort.list(
                myCodonPositions$position,
                decreasing = FALSE
            ), ]
        rownames(myCodonPositions) <- NULL

        return(myCodonPositions)
    }
    codonsOfInterest <-
        myFindPotentialStartsAndStops(
            dnaSequence = dnaSequence,
            codons = codonAnnotation$codons
        )

    ### Outcommented 19/3/21 to discard identification of truncated ORFs
    ## Add stop codon at the end to make sure ORFs are allowed to continue over the edge of the transcript (by simmulating the last codons in each reading frame is a stop codon)
    #codonsOfInterest <-
    #    rbind(
    #        codonsOfInterest,
    #        data.frame(
    #            position = nchar(dnaSequence) - 2 - 2:0,
    #            codon = codonAnnotation$codons[which(
    #                codonAnnotation$meaning == 'stop'
    #            )][1],
    #            stringsAsFactors = FALSE
    #        )
    #    )

    ### Annotate with meaning of condon
    codonsOfInterest$meaning <-
        codonAnnotation$meaning[match(
            codonsOfInterest$codon,
            codonAnnotation$codons
        )]

    # Filter for annotated start sites
    if (!is.null(filterForPostitions)) {
        if( any( !is.na( filterForPostitions) ) ) { # KVS added Dec 2020 : NAs are how those without overlapping CDS was added back for "longest.AnnotatedWhenPossible"
            codonsOfInterest <-
                codonsOfInterest[which(
                    codonsOfInterest$position %in% filterForPostitions |
                        codonsOfInterest$meaning != 'start'
                ), ]
            if (!nrow(codonsOfInterest)) {
                return(data.frame(NULL))
            }
        }
    }


    # Loop over the 3 possible reading frames
    myORFs <- list()
    for (i in 0:2) {
        # Reduce the data.frame with codons to only contain those within that reading frame
        localCodonsOfInterest <-
            codonsOfInterest[which(codonsOfInterest$position %% 3 == i), ]

        # if there are any look for ORFs
        if (nrow(localCodonsOfInterest) == 0) {
            next
        }

        ### Create rle objet to allow analysis of order for start/stop codons
        myRle <- rle(localCodonsOfInterest$meaning)

        ### Trim codons
        # Remove rows so the first codon is a start codon (not a stop codon)
        redoRLE <- FALSE
        if (myRle$values[1] == 'stop') {
            localCodonsOfInterest <-
                localCodonsOfInterest[
                    (myRle$lengths[1] + 1):(nrow(localCodonsOfInterest))
                , ]

            #redo rle nessesary
            redoRLE <- TRUE
        }
        # Remove rows so the last codon is a stop codon (not a start codon)
        if (tail(myRle$values, 1) == 'start') {
            localCodonsOfInterest <-
                localCodonsOfInterest[
                    1:(nrow(localCodonsOfInterest) - tail(myRle$lengths, 1))
                , ]

            #redo rle nessesary
            redoRLE <- TRUE
        }
        # redo RLE if anything was removed
        if (redoRLE) {
            myRle <- rle(localCodonsOfInterest$meaning)
        }

        # make sure that there are both start and stop (left)
        if (!all(c('start', 'stop') %in% localCodonsOfInterest$meaning)) {
            next
        }

        # Remove duplicates to make sure that I only take the first start codon (if multiple are pressent) and the first stop codon if multiple are pressent
        localCodonsOfInterest <-
            localCodonsOfInterest[
                c(1, cumsum(myRle$lengths) + 1)[1:length(myRle$lengths)]
            , ]

        ### Loop over the resulting table and concattenate the ORFs
        localCodonsOfInterest$myOrf <-
            as.vector(sapply(1:(nrow(
                localCodonsOfInterest
            ) / 2), function(x)
                rep(x, 2)))
        localCodonsOfInterestSplit <- split(
                localCodonsOfInterest$position,
                f = localCodonsOfInterest$myOrf
        )
        myORFs[[i + 1]] <-
            list(myListToDf(
                lapply(localCodonsOfInterestSplit, function(x)
                    data.frame(start = x[1], end = x[2]))
            ))
    }
    if (length(myORFs) == 0) {
        return(data.frame(NULL))
    }

    # Massage from list to data.frame
    myORFs <-
        myORFs[which(!sapply(myORFs, is.null))] # remove potential empty ones
    myORFs <-
        myListToDf(lapply(myORFs, function(x)
            x[[1]])) # x[[1]] nessary because of the extra list i had to introduce to avoid classes due to single entry lists

    ### Make final calculatons
    # Subtract 1 since the positions I here have worked with are the start of the stop codon positions (and I do not want the stop codon to be included)
    myORFs$end <- myORFs$end - 1
    # add lengths
    myORFs$length <-
        myORFs$end - myORFs$start + 1 # the +1 is because both are included

    # sort
    myORFs <- myORFs[sort.list(myORFs$start, decreasing = FALSE), ]
    rownames(myORFs) <- NULL

    return(myORFs)
}

analyzeORFforPTC <- function(
    aORF,           # A data.frame containing the transcript coordinats and length of the main ORF
    exonStructure,  # A data.frame with the exon structure of the transcript
    PTCDistance = 50  # A numeric giving the premature termination codon-distance: The minimum distance from a STOP to the final exon-exon junction, for a transcript to be marked as NMD-sensitive
){
    # Have to be done diffetly for each strand due to the reversed exon structure
    if (exonStructure$strand[1] == '+') {
        # Calculate exon cumSums (because they "translate" the genomic coordinats to transcript coordinats )
        exonCumsum      <- cumsum(c(0,      exonStructure$width))

        # Calculate wich exon the start and stop codons are in
        cdsStartExonIndex   <- max(which(aORF$start >  exonCumsum))
        cdsEndExonIndex     <- max(which(aORF$end   >  exonCumsum))
        # Calcualte genomic position of the ORF
        cdsStartGenomic <-
            exonStructure$start[cdsStartExonIndex]  +
            (aORF$start - exonCumsum[cdsStartExonIndex]  - 1) # -1 because both intervals are inclusive [a;b] (meaning the start postition will be counted twice)
        cdsEndGenomic   <-
            exonStructure$start[cdsEndExonIndex]    +
            (aORF$end   - exonCumsum[cdsEndExonIndex]  - 1) # -1 because both intervals are inclusive [a;b] (meaning the start postition will be counted twice)

        # Calculate length from stop codon to last splice junction
        #               transcript pos for last exon:exon junction  - transcript pos for stop codon
        stopDistance <-
            exonCumsum[length(exonCumsum) - 1]  -
            aORF$end # (positive numbers means the stop position upstream of the last exon exon junction )

        # Calculate which exon the stop codon are in compared to the last exon exon junction
        # stop in exon      total nr exon
        junctionIndex <- cdsEndExonIndex - nrow(exonStructure)
    }
    if (exonStructure$strand[1] == '-') {
        # Calculate exon cumSums (because they "translate" the genomic coordinats to transcript coordinats )
        exonRevCumsum   <- cumsum(c(0, rev(exonStructure$width)))

        # Calculate wich exon the start and stop codons are in
        cdsStartExonIndex   <-
            max(which(aORF$start >  exonRevCumsum))
        cdsEndExonIndex     <-
            max(which(aORF$end   >  exonRevCumsum))

        # create a vector to translate indexes to reverse (needed when exon coordinats are extracted)
        reversIndexes <- nrow(exonStructure):1

        # Calcualte genomic position of the ORF (end and start are switched in order to return them so start < end (default of all formating))
        cdsStartGenomic <-
            exonStructure$end[reversIndexes[cdsStartExonIndex]]  -
            (aORF$start - exonRevCumsum[cdsStartExonIndex] - 1) # -1 because both intervals are inclusive [a;b] (meaning the start postition will be counted twice)
        cdsEndGenomic   <-
            exonStructure$end[reversIndexes[cdsEndExonIndex]]  -
            (aORF$end   - exonRevCumsum[cdsEndExonIndex] - 1) # -1 because both intervals are inclusive [a;b] (meaning the start postition will be counted twice)

        # Calculate length from stop codon to last splice junction
        #               transcript pos for last exon:exon junction  - transcript pos for stop codon
        stopDistance <-
            exonRevCumsum[length(exonRevCumsum) - 1]    - aORF$end # (positive numbers means the stop position upstream of the last exon exon junction )

        # Calculate which exon the stop codon are in compared to the last exon exon junction
        # stop in exon      total nr exon
        junctionIndex <- cdsEndExonIndex - nrow(exonStructure)
    }

    return(
        data.frame(
            orfStarExon = cdsStartExonIndex,
            orfEndExon = cdsEndExonIndex,
            orfStartGenomic = cdsStartGenomic,
            orfEndGenomic = cdsEndGenomic,
            stopDistanceToLastJunction = stopDistance,
            stopIndex = junctionIndex,
            PTC = (stopDistance >= PTCDistance) & junctionIndex != 0
        )
    )

}
