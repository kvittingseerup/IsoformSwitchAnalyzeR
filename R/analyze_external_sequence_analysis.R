##############################
### Import and add external analysis

# Helper function
convertCoordinatsTranscriptToGenomic <- function(
    transcriptCoordinats,   # A data.frame containing the transcript coordinats
    exonStructure           # A data.frame with the genomic coordinats of the transcript exons
) {
    if (exonStructure$strand[1] == '+') {
        # Calculate exon cumSums (because they "translate" the genomic coordinats to transcript coordinats )
        exonCumsum      <- cumsum(c(0,      exonStructure$width))

        # Calculate wich exon the start and stop codons are in
        cdsStartExonIndex   <-
            max(which(transcriptCoordinats$start >  exonCumsum))
        cdsEndExonIndex     <-
            max(which(transcriptCoordinats$end   >  exonCumsum))
        # Calcualte genomic position of the ORF
        cdsStartGenomic <-
            exonStructure$start[cdsStartExonIndex]  +
            (transcriptCoordinats$start - exonCumsum[cdsStartExonIndex]  - 1) # -1 because both intervals are inclusive [a;b] (meaning the start postition will be counted twice)
        cdsEndGenomic   <-
            exonStructure$start[cdsEndExonIndex]    +
            (transcriptCoordinats$end   - exonCumsum[cdsEndExonIndex]  - 1) # -1 because both intervals are inclusive [a;b] (meaning the start postition will be counted twice)

    }
    if (exonStructure$strand[1] == '-') {
        # Calculate exon cumSums (because they "translate" the genomic coordinats to transcript coordinats )
        exonRevCumsum   <- cumsum(c(0, rev(exonStructure$width)))

        # Calculate wich exon the start and stop codons are in
        cdsStartExonIndex   <-
            max(which(transcriptCoordinats$start >  exonRevCumsum))
        cdsEndExonIndex     <-
            max(which(transcriptCoordinats$end   >  exonRevCumsum))

        # create a vector to translate indexes to reverse (needed when exon coordinats are extracted)
        reversIndexes <- nrow(exonStructure):1

        # Calcualte genomic position of the ORF (end and start are switched in order to return them so start < end (default of all formating))
        cdsStartGenomic <-
            exonStructure$end[reversIndexes[cdsStartExonIndex]]  -
            (transcriptCoordinats$start - exonRevCumsum[cdsStartExonIndex] - 1) # -1 because both intervals are inclusive [a;b] (meaning the start postition will be counted twice)
        cdsEndGenomic   <-
            exonStructure$end[reversIndexes[cdsEndExonIndex]]  -
            (transcriptCoordinats$end   - exonRevCumsum[cdsEndExonIndex] - 1) # -1 because both intervals are inclusive [a;b] (meaning the start postition will be counted twice)
    }

    return(
        data.frame(
            pfamStarExon = cdsStartExonIndex,
            pfamEndExon = cdsEndExonIndex,
            pfamStartGenomic = cdsStartGenomic,
            pfamEndGenomic = cdsEndGenomic
        )
    )

}


### Actural functions

analyzeCPAT <- function(
    switchAnalyzeRlist,
    pathToCPATresultFile,
    codingCutoff,
    removeNoncodinORFs,
    quiet = FALSE
) {
    ### Check input
    if (TRUE) {
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(
                'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
            )
        }

        # file
        if (class(pathToCPATresultFile) != 'character') {
            stop(
                'The \'pathToCPATresultFile\' argument must be a string pointing to the CPAT result file'
            )
        }
        if (!file.exists(pathToCPATresultFile)) {
            stop('The file \'pathToCPATresultFile\' points to does not exist')
        }

        # codingCutoff
        if (missing(codingCutoff)) {
            stop('The \'codingCutoff\' argument must be supplied')
        }
        if (length(codingCutoff) != 1 |
            class(codingCutoff) != 'numeric') {
            stop('The \'codingCutoff\' argument must be a numeric of length 1 ')
        }
        if (!codingCutoff >= 0 &
            codingCutoff <= 1) {
            stop('The \'codingCutoff\' argument must be a numeric in the inverval [0,1]')
        }

        if (missing(removeNoncodinORFs)) {
            stop('The \'removeNoncodinORFs\' argument must be supplied')
        }
        if (!is.logical(removeNoncodinORFs)) {
            stop('The \'removeNoncodinORFs\' argument must be a logical (TRUE or FALSE)')
        }
    }

    ### Read file
    if (TRUE) {
        ### read in file
        myCPATresults <-
            read.table(
                file = pathToCPATresultFile,
                stringsAsFactors = FALSE,
                header = TRUE,
                sep = '\t'
            )

        if( nrow(myCPATresults) == 0) {
            stop('No results were found in the result file')
        }

        # check if it is web file
        if (ncol(myCPATresults) == 8) {
            temp <- myCPATresults$Sequence.Name
            myCPATresults <- myCPATresults[, 3:7]
            myCPATresults$id <- temp

            if (!all(
                colnames(myCPATresults) %in% c(
                    "RNA.size",
                    "ORF.size",
                    "Ficket.Score",
                    "Hexamer.Score",
                    "Coding.Probability",
                    "id"
                )
            )) {
                stop(
                    'There seems to be a problem with the CPAT result file. Please check it is the rigth file and try again'
                )
            }
            # rename to match non-web file
            colnames(myCPATresults) <-
                c(
                    'mRNA_size',
                    'ORF_size',
                    'Fickett_score',
                    'Hexamer_score',
                    'coding_prob',
                    'id'
                )


        } else if (ncol(myCPATresults) == 5) {
            myCPATresults$id <- rownames(myCPATresults)
            rownames(myCPATresults) <- NULL

            # Check file header
            if (!all(
                colnames(myCPATresults) == c(
                    'mRNA_size',
                    'ORF_size',
                    'Fickett_score',
                    'Hexamer_score',
                    'coding_prob',
                    'id'
                )
            )) {
                stop(
                    'The file pointed to by the \'pathToCPATresultFile\' argument is not a CPAT result file'
                )
            }


        } else {
            stop(
                'There seems to be a problem with the CPAT result file. Please check it is the rigth file and try again'
            )
        }

        ### Massage
        # check ids
        if (!any(
            tolower(myCPATresults$id) %in%
            tolower(switchAnalyzeRlist$isoformFeatures$isoform_id)
        )) {
            stop(
                'The transcript ids in the file pointed to by the \'pathToCPATresultFile\' argument does not match the transcripts stored in the supplied switchAnalyzeRlist'
            )
        }

        ### Overwirte with correct cased names
        myCPATresults$id <-
            switchAnalyzeRlist$isoformFeatures$isoform_id[match(
                tolower(myCPATresults$id),
                tolower(switchAnalyzeRlist$isoformFeatures$isoform_id)
            )]

        ### subset
        myCPATresults <-
            myCPATresults[which(
                myCPATresults$id %in%
                    switchAnalyzeRlist$isoformFeatures$isoform_id
            ), ]
    }

    ### Add analysis to switchAnalyzeRlist
    if (TRUE) {
        matchVector <-
            match(switchAnalyzeRlist$isoformFeatures$isoform_id,
                  myCPATresults$id)

        switchAnalyzeRlist$isoformFeatures$codingPotentialValue <-
            myCPATresults$coding_prob[matchVector]
        switchAnalyzeRlist$isoformFeatures$codingPotential      <-
            myCPATresults$coding_prob[matchVector] >= codingCutoff


        n <- length(unique(myCPATresults$id))
        p <-
            round(n / length(
                unique(switchAnalyzeRlist$isoformFeatures$isoform_id)
            ) * 100, digits = 2)

        if (!quiet) {
            message(paste(
                'Added coding potential to ',
                n,
                ' (',
                p,
                '%) transcripts',
                sep = ''
            ))
        }
    }

    if (removeNoncodinORFs &
        !is.null(switchAnalyzeRlist$orfAnalysis)) {
        ### Extract coding isoforms
        nonCodingIsoforms <-
            unique(switchAnalyzeRlist$isoformFeatures$isoform_id[which(
                !switchAnalyzeRlist$isoformFeatures$codingPotential
            )])

        ### Replace noncoding isoforms
        switchAnalyzeRlist$orfAnalysis[which(
            switchAnalyzeRlist$orfAnalysis$isoform_id %in% nonCodingIsoforms),
            2:ncol(switchAnalyzeRlist$orfAnalysis)
        ] <- NA

        ### Overwrite PTC
        switchAnalyzeRlist$isoformFeatures$PTC[which(
            switchAnalyzeRlist$isoformFeatures$isoform_id %in%
                nonCodingIsoforms
        )] <- NA
    }

    return(switchAnalyzeRlist)
}


analyzePFAM <- function(
    switchAnalyzeRlist,
    pathToPFAMresultFile,
    showProgress = TRUE,
    quiet = FALSE
) {
    if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
        stop(
            'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
        )
    }
    if (is.null(switchAnalyzeRlist$orfAnalysis)) {
        stop('ORF needs to be analyzed. Please run analyzeORF and try again.')
    }

    # file
    if (class(pathToPFAMresultFile) != 'character') {
        stop(
            'The \'pathToPFAMresultFile\' argument must be a string pointing to the PFAM result file'
        )
    }
    if (!file.exists(pathToPFAMresultFile)) {
        stop('The file \'pathToPFAMresultFile\' points to does not exist')
    }


    ### Test wither headers are included
    temp <-
        read.table(
            file = pathToPFAMresultFile,
            stringsAsFactors = FALSE,
            fill = TRUE,
            header = FALSE,
            nrows = 1
        )
    if (nrow(temp) == 0) {
        stop('The file pointed to by \'pathToPFAMresultFile\' is empty')
    }
    if (grepl('^<seq|^seq', temp[1, 1])) {
        skipLine <- 1
    } else {
        skipLine <- 0
    }

    ### read in pfam resut result
    myPfamResult <-
        read.table(
            file = pathToPFAMresultFile,
            stringsAsFactors = FALSE,
            fill = TRUE,
            header = FALSE,
            col.names = 1:16,
            skip = skipLine
        )
    if (nrow(myPfamResult) == 0) {
        stop('The file pointed to by \'pathToPFAMresultFile\' is empty')
    }

    ### Identify which type of PFAM file is supplied
    # Colnames for pfam result
    oldColnames <-
        c(
            'seq_id',
            'alignment_start',
            'alignment_end',
            'envelope_start',
            'envelope_end',
            'hmm_acc',
            'hmm_name',
            'type',
            'hmm_start',
            'hmm_end',
            'hmm_length',
            'bit_score',
            'E_value',
            'significant',
            'clan',
            'residue'
        )
    newColNames <-
        c(
            'seq_id',
            'alignment_start',
            'alignment_end',
            'envelope_start',
            'envelope_end',
            'hmm_acc',
            'hmm_name',
            'hmm_start',
            'hmm_end',
            'hmm_length',
            'bit_score',
            'Individual_E_value',
            'Conditional_E_value',
            'significant',
            'outcompeted',
            'clan'
        )

    ### Old style
    if (class(myPfamResult$X8) == 'character') {
        colnames(myPfamResult) <- oldColnames
        myPfamResult$residue[which(myPfamResult$residue == '')] <-
            NA

        ### New style
    } else if (class(myPfamResult$X8) == 'integer') {
        colnames(myPfamResult) <- newColNames
        myPfamResult$clan[which(myPfamResult$clan == '')] <- NA

    } else {
        stop('The file supplied is not recogniced as a pfam output.')
    }

    ### Sanity check that it is a PFAM result file
    if (TRUE) {
        test1 <-
            ncol(myPfamResult) == 15 |
            ncol(myPfamResult) == 16 # the output have 15 or 16 collumns depending on whther active sites are predicted
        test2 <-
            all(grepl(
                pattern = '^PF|^PB' ,
                myPfamResult$hmm_acc,
                ignore.case = FALSE
            ))                # All pfam hmm starts with PF

        if (!all(test1, test2)) {
            stop('The file pointed to by \'pathToPFAMresultFile\' is not a PFAM result file')
        }

        # test names
        if (!any(myPfamResult$seq_id %in%
                 switchAnalyzeRlist$isoformFeatures$isoform_id)) {
            stop(
                'The transcript ids in the file pointed to by the \'pathToPFAMresultFile\' argument does not match the transcripts stored in the supplied switchAnalyzeRlist'
            )
        }

        ### Only add to those with annotated ORF
        #myPfamResult <- myPfamResult[which( myPfamResult$seq_id %in% switchAnalyzeRlist$isoformFeatures$isoform_id ),]
        myPfamResult <- myPfamResult[which(
            myPfamResult$seq_id %in%
                switchAnalyzeRlist$orfAnalysis$isoform_id[which(
                    !is.na(switchAnalyzeRlist$orfAnalysis$orfTransciptStart)
                )]
        ), ]
    }

    if (showProgress & !quiet) {
        progressBar <- 'text'
    } else {
        progressBar <- 'none'
    }

    ### Fill in blanks if active residues are included
    if (TRUE) {
        if ('residue' %in% colnames(myPfamResult)) {
            # test whether active residue analysis was performed
            if (any(nchar(myPfamResult$residue) != 0) &
                all(!is.na(myPfamResult$residue))) {
                withActiveRes <- TRUE
            } else {
                withActiveRes <- FALSE
            }

            # If active residue analysis was performed do the additional massage of the data
            if (withActiveRes) {
                ### Exchange empty strings in residus with NA's
                myPfamResult$residue[which(!grepl(
                    'predicted_active_site',
                    myPfamResult$residue
                ))] <- NA
                colnames(myPfamResult)[which(
                    colnames(myPfamResult) == 'residue')] <-
                    'predicted_active_site'

                ### convert active residue to somthing readable
                activeSiteIndex <-
                    which(!is.na(myPfamResult$predicted_active_site))
                myPfamResult$predicted_active_site[activeSiteIndex] <-
                    sapply(
                        myPfamResult$predicted_active_site[activeSiteIndex],
                        function(aVec) {
                            substr(aVec,
                                   start = 23,
                                   stop = nchar(aVec) - 1)
                        }
                    )
            }
        } else {
            withActiveRes <- FALSE
        }


    }

    ### Convert from AA coordinats to transcript and genomic coordinats
    if (TRUE) {
        if (!quiet) {
            message('Converting AA coordinats to transcript and genomic coordinats...')
        }
        ### Remove unwanted columns
        myPfamResult$envelope_start <- NULL
        myPfamResult$envelope_end <- NULL
        myPfamResult$hmm_start <- NULL
        myPfamResult$hmm_end <- NULL
        myPfamResult$hmm_length <- NULL

        colnames(myPfamResult)[which(
            grepl('alignment_', colnames(myPfamResult))
        )] <- c('orf_aa_start', 'orf_aa_end')

        ### convert from codons to transcript position
        orfStartDF <-
            unique(as.data.frame(
                switchAnalyzeRlist$orfAnalysis[,
                    c('isoform_id', 'orfTransciptStart')
                ]
            ))
        myPfamResult$transcriptStart <-
            (myPfamResult$orf_aa_start  * 3 - 2) +
            orfStartDF[
                match(
                    x = myPfamResult$seq_id,
                    table = orfStartDF$isoform_id
                ),
                2] - 1
        myPfamResult$transcriptEnd <-
            (myPfamResult$orf_aa_end * 3) +
            orfStartDF[
                match(
                    x = myPfamResult$seq_id,
                    table = orfStartDF$isoform_id
                ),
                2] - 1

        ### convert from transcript to genomic coordinats
        # extract exon data
        myExons <-
            as.data.frame(switchAnalyzeRlist$exons[which(
                switchAnalyzeRlist$exons$isoform_id %in% myPfamResult$seq_id
            ), ])
        myExonsSplit <- split(myExons, f = myExons$isoform_id)

        # loop over the individual transcripts and extract the genomic coordiants of the domain and also for the active residues (takes 2 min for 17000 rows)
        myPfamResultDf <-
            plyr::ddply(
                myPfamResult,
                .progress = progressBar,
                .variables = 'seq_id',
                .fun = function(aDF) {
                    transcriptId <- aDF$seq_id[1]
                    localExons <-
                        as.data.frame(myExonsSplit[[transcriptId]])

                    # extract domain allignement
                    localORFalignment <- aDF
                    colnames(localORFalignment)[match(
                        x = c('transcriptStart', 'transcriptEnd'),
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

                        ### look into active residues
                        if (withActiveRes) {
                            if (!is.na(
                                localORFalignment$predicted_active_site[j]
                            )) {
                                activeResInfo <-
                                    data.frame(activeRes = as.integer(unlist(
                                        strsplit(
                                            x = localORFalignment$predicted_active_site[j],
                                            split = ','
                                        )
                                    )))

                                activeResInfo$start <-
                                    activeResInfo$activeRes  * 3 - 2
                                activeResInfo$end   <-
                                    activeResInfo$activeRes  * 3

                                activeResInfoList <- list()
                                for (k in 1:nrow(activeResInfo)) {
                                    activeResInfoList[[as.character(k)]] <-
                                        convertCoordinatsTranscriptToGenomic(
                                            transcriptCoordinats = activeResInfo[k, ],
                                            exonStructure = localExons
                                        )[, c('pfamStartGenomic',
                                              'pfamEndGenomic')]
                                }
                                activeResInfoDf <-
                                    cbind(activeResInfo,
                                          do.call(rbind, activeResInfoList))

                                ### add it to the domain info
                                domainInfo$activeResTranscriptStart <-
                                    paste(activeResInfoDf$start, collapse = ',')
                                domainInfo$activeResTranscriptEnd   <-
                                    paste(activeResInfoDf$end, collapse = ',')
                                domainInfo$activeResGenomicStart    <-
                                    paste(activeResInfoDf$pfamStartGenomic,
                                          collapse = ',')
                                domainInfo$activeResGenomicEnd      <-
                                    paste(activeResInfoDf$pfamEndGenomic,
                                          collapse = ',')
                            } else {
                                ### Add NA instead of residues
                                domainInfo$activeResTranscriptStart <- NA
                                domainInfo$activeResTranscriptEnd   <- NA
                                domainInfo$activeResGenomicStart    <- NA
                                domainInfo$activeResGenomicEnd      <- NA
                            }
                        }


                        orfPosList[[as.character(j)]] <- domainInfo

                    }
                    orfPosDf <- do.call(rbind, orfPosList)

                    return(cbind(aDF, orfPosDf))
                }
            )

    }

    ### Add analysis to switchAnalyzeRlist
    if (TRUE) {
        ### reorder data.frame
        # make sure the basic data is last
        colnames(myPfamResultDf)[1] <- 'isoform_id'
        newOrderNames <- c('isoform_id', 'hmm_acc', 'hmm_name', 'clan')
        myPfamResultDf <-
            myPfamResultDf[, c(
                which(colnames(myPfamResultDf) %in% newOrderNames) ,
                which(!colnames(myPfamResultDf) %in% newOrderNames)
            )]

        # sort
        myPfamResultDf <-
            myPfamResultDf[order(
                myPfamResultDf$isoform_id,
                myPfamResultDf$transcriptStart,
                myPfamResultDf$hmm_name
            ), ]

        # if active residues are pressent put them last
        if (withActiveRes) {
            newOrder <-
                c(which(
                    !grepl(
                        'Predicted_active_site|ActiveRes',
                        colnames(myPfamResultDf)
                    )
                ), which(
                    grepl(
                        'Predicted_active_site|ActiveRes',
                        colnames(myPfamResultDf)
                    )
                ))
            myPfamResultDf <- myPfamResultDf[, newOrder]
        }

        # add the pfam results to the switchAnalyzeRlist object
        switchAnalyzeRlist$domainAnalysis <- myPfamResultDf

        # add indication to transcriptDf
        switchAnalyzeRlist$isoformFeatures$domain_identified <- 'no'
        switchAnalyzeRlist$isoformFeatures$domain_identified[which(
            is.na(switchAnalyzeRlist$isoformFeatures$PTC)
        )] <- NA # sets NA for those not analyzed

        switchAnalyzeRlist$isoformFeatures$domain_identified[which(
            switchAnalyzeRlist$isoformFeatures$isoform_id %in%
                myPfamResultDf$isoform_id
        )] <- 'yes'
    }

    n <- length(unique(myPfamResultDf$isoform_id))
    p <-
        round(n / length(unique(
            switchAnalyzeRlist$isoformFeatures$isoform_id
        )) * 100, digits = 2)

    if (!quiet) {
        message(paste(
            'Added domain information to ',
            n,
            ' (',
            p,
            '%) transcripts',
            sep = ''
        ))
    }
    return(switchAnalyzeRlist)
}


analyzeSignalP <- function(
    switchAnalyzeRlist,
    pathToSignalPresultFile,
    quiet = FALSE
) {
    if (is.null(switchAnalyzeRlist$orfAnalysis)) {
        stop('ORF needs to be analyzed. Please run analyzeORF and try again.')
    }

    # file
    if (class(pathToSignalPresultFile) != 'character') {
        stop(
            'The \'pathToSignalPresultFile\' argument must be a string pointing to the SignalP result file'
        )
    }
    if (!file.exists(pathToSignalPresultFile)) {
        stop('The file \'pathToSignalPresultFile\' points to does not exist')
    }


    ### Obtain signalP result
    if (TRUE) {
        ### Read file
        singalPresults <-
            read.table(
                pathToSignalPresultFile,
                header = FALSE,
                stringsAsFactors = FALSE,
                fill = TRUE,
                col.names = paste('V', 1:13, sep = '')
            )
        if(nrow(singalPresults) == 0) {
            stop('The result file seems to be empty')
        }

        # extract summary
        singalPresults <-
            singalPresults[which(grepl(
                pattern = "SP=\'YES\'",
                x =  singalPresults$V2
            )), ]

        ### Sanity check that it is a SignalIP result file
        if (TRUE) {
            if (nrow(singalPresults) == 0) {
                stop('No signial peptides were found')
            }

            t1 <- all(singalPresults$V3 == 'Cleavage')
            t2 <- all(singalPresults$V4 == 'site')
            t3 <- all(singalPresults$V5 == 'between')
            t4 <- all(singalPresults$V6 == 'pos.')
            t5 <- all(grepl('^Name=', singalPresults$V1))
            t6 <- all(grepl('^Networks=', singalPresults$V13))

            if (!all(t1, t2, t3, t4, t5, t6)) {
                stop(
                    'The file pointed to by \'pathToSignalPresultFile\' is not a SignalIP summary file'
                )
            }
        }
    }

    ### Convert from AA coordinats to transcript and genomic coordinats
    if (TRUE) {
        ### Massage data
        singalPresults <- singalPresults[, c(1, 7, 13)]
        colnames(singalPresults) <-
            c('transcript_id', 'aa_removed', 'network_used')
        singalPresults$transcript_id <-
            sapply(
                strsplit(singalPresults$transcript_id, 'Name='),
                function(x) {x[2]}
            )
        singalPresults$network_used <-
            sapply(
                strsplit(singalPresults$network_used, 'Networks=SignalP-'),
                function(x) {
                    x[2]
                }
            )

        # test names
        if (!any(
            singalPresults$transcript_id %in%
            switchAnalyzeRlist$isoformFeatures$isoform_id
        )) {
            stop(
                'The transcript ids in the file pointed to by the \'pathToSignalPresultFile\' argument does not match the transcripts stored in the supplied switchAnalyzeRlist'
            )
        }

        ### Reduce to relevant
        singalPresults <- singalPresults[which(
            singalPresults$transcript_id %in%
                switchAnalyzeRlist$orfAnalysis$isoform_id[which(
                    !is.na(switchAnalyzeRlist$orfAnalysis$orfTransciptStart)
                )]
        ), ]

        ### Convert cleavage site to transcript coordinats
        # df with orf start
        orfStartDF <-
            unique(as.data.frame(switchAnalyzeRlist$orfAnalysis[,c(
                'isoform_id', 'orfTransciptStart'
            )]))
        # calculate start position
        singalPresults$transcriptClevageAfter <-
            (singalPresults$aa_removed  * 3) + orfStartDF[match(
                x = singalPresults$transcript_id,
                table = orfStartDF$isoform_id
            ), 2] - 1

        ### convert from transcript to genomic coordinats
        # extract exon data
        myExons <-
            as.data.frame(switchAnalyzeRlist$exons[which(
                switchAnalyzeRlist$exons$isoform_id %in%
                    singalPresults$transcript_id
            ),])
        myExonsSplit <- split(myExons, f = myExons$isoform_id)

        # calculate genomic coordinat
        singalPresults <-
            plyr::ddply(
                singalPresults,
                .variables = 'transcript_id',
                .fun = function(aDF) {
                    localExons <- myExonsSplit[[aDF$transcript_id[1]]]

                    tempDf <-
                        data.frame(
                            start = aDF$transcriptClevageAfter,
                            end = aDF$transcriptClevageAfter
                        )
                    aDF$genomicClevageAfter <-
                        convertCoordinatsTranscriptToGenomic(
                            transcriptCoordinats = tempDf,
                            exonStructure = localExons
                        )$pfamStartGenomic
                    return(aDF)
                }
            )

        singalPresults$has_signal_peptide <- 'yes'
        singalPresults <-
            singalPresults[, c(
                'transcript_id',
                'has_signal_peptide',
                'network_used',
                'aa_removed',
                'transcriptClevageAfter',
                'genomicClevageAfter'
            )]

        colnames(singalPresults) <-
            gsub('transcript_id',
                 'isoform_id',
                 colnames(singalPresults))

    }

    ### Add analysis to switchAnalyzeRlist
    if (TRUE) {
        # add the pfam results to the switchAnalyzeRlist object
        switchAnalyzeRlist$signalPeptideAnalysis <- singalPresults

        # add indication to transcriptDf
        switchAnalyzeRlist$isoformFeatures$signal_peptide_identified <-
            'no'
        switchAnalyzeRlist$isoformFeatures$signal_peptide_identified[which(
            is.na(switchAnalyzeRlist$isoformFeatures$PTC)
        )] <- NA # sets NA for those not analyzed
        switchAnalyzeRlist$isoformFeatures$signal_peptide_identified[which(
            switchAnalyzeRlist$isoformFeatures$isoform_id %in%
                singalPresults$isoform_id
        )] <- 'yes'

        n <- length(unique(singalPresults$isoform_id))
        p <-
            round(n / length(
                unique(switchAnalyzeRlist$isoformFeatures$isoform_id)
            ) * 100, digits = 2)
        if (!quiet) {
            message(
                paste(
                    'Added signal peptide information to ',
                    n,
                    ' (',
                    p,
                    '%) transcripts',
                    sep = ''
                )
            )
        }
    }

    return(switchAnalyzeRlist)
}
