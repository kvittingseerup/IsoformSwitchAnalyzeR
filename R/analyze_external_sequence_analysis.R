### Run external analysis

##############################
### Import and add external analysis

### Actural functions
analyzeCPAT <- function(
    ### Core arguments
    switchAnalyzeRlist,
    pathToCPATresultFile,
    codingCutoff,
    removeNoncodinORFs,
    
    ### Advanced arguments
    ignoreAfterBar = TRUE,
    ignoreAfterSpace = TRUE,
    ignoreAfterPeriod = FALSE,
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
            stop('No results were found in the CPAT result file')
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
        
        ### Fix names
        if( ignoreAfterBar | ignoreAfterSpace | ignoreAfterPeriod){
          myCPATresults$id <- fixNames(
            nameVec = myCPATresults$id,
            ignoreAfterBar = ignoreAfterBar,
            ignoreAfterSpace = ignoreAfterSpace,
            ignoreAfterPeriod = ignoreAfterPeriod
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
        toModify <- which( switchAnalyzeRlist$orfAnalysis$isoform_id %in% nonCodingIsoforms)
        switchAnalyzeRlist$orfAnalysis[
            toModify,
            which( ! colnames(switchAnalyzeRlist$orfAnalysis) %in% c('isoform_id','orf_origin'))
        ] <- NA
        if( ! is.null(switchAnalyzeRlist$orfAnalysis$orf_origin)) {
            switchAnalyzeRlist$orfAnalysis$orf_origin[toModify] <- 'Predicted'
        }

        ### Overwrite PTC
        switchAnalyzeRlist$isoformFeatures$PTC[which(
            switchAnalyzeRlist$isoformFeatures$isoform_id %in%
                nonCodingIsoforms
        )] <- NA
    }

    return(switchAnalyzeRlist)
}

analyzeCPC2 <- function(
    ### Core arguments
    switchAnalyzeRlist,
    pathToCPC2resultFile,
    codingCutoff = 0.5,
    removeNoncodinORFs,
    
    ### Advanced arguments
    ignoreAfterBar = TRUE,
    ignoreAfterSpace = TRUE,
    ignoreAfterPeriod = FALSE,
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
        if (class(pathToCPC2resultFile) != 'character') {
            stop(
                'The \'pathToCPC2resultFile\' argument must be a string pointing to the CPC2 result file'
            )
        }
        if (!file.exists(pathToCPC2resultFile)) {
            stop('The file \'pathToCPC2resultFile\' points to does not exist')
        }

        # codingCutoff
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
        myCPCresults <-
            read.table(
                file = pathToCPC2resultFile,
                stringsAsFactors = FALSE,
                header = TRUE,
                comment.char = '',
                sep = '\t'
            )

        if( nrow(myCPCresults) == 0) {
            stop('No results were found in the CPC result file')
        }

        # check if it is web file
        if (ncol(myCPCresults) == 7) {
            colsToMatch <- c(
                "X.ID",
                "peptide_length",
                "Fickett_score",
                "pI",
                "ORF_integrity",
                "coding_probability",
                'label'
            )

            if (!all(
                colsToMatch == colnames(myCPCresults)
            )) {
                stop(
                    'There seems to be a problem with the CPC2 result file. Please check it is the rigth file and try again'
                )
            }

        } else if (ncol(myCPCresults) == 8) {
            colsToMatch <- c(
                "X.ID",
                'transcript_length',
                "peptide_length",
                "Fickett_score",
                "pI",
                "ORF_integrity",
                "coding_probability",
                'label'
            )

            if (!all(
                colsToMatch == colnames(myCPCresults)
            )) {
                stop(
                    'There seems to be a problem with the CPC2 result file. Please check it is the rigth file and try again'
                )
            }

        } else {
            stop(
                'There seems to be a problem with the CPC2 result file. Please check it is the rigth file and try again'
            )
        }
        # rename to match non-web file
        colnames(myCPCresults)[1] <- 'id'
        
        ### Fix names
        if( ignoreAfterBar | ignoreAfterSpace | ignoreAfterPeriod){
          myCPCresults$id <- fixNames(
            nameVec = myCPCresults$id,
            ignoreAfterBar = ignoreAfterBar,
            ignoreAfterSpace = ignoreAfterSpace,
            ignoreAfterPeriod = ignoreAfterPeriod
          )
        }

        ### Massage
        # check ids
        if (!any(
            tolower(myCPCresults$id) %in%
            tolower(switchAnalyzeRlist$isoformFeatures$isoform_id)
        )) {
            stop(
                'The transcript ids in the file pointed to by the \'pathToCPC2resultFile\' argument does not match the transcripts stored in the supplied switchAnalyzeRlist'
            )
        }

        ### Overwirte with correct cased names
        myCPCresults$id <-
            switchAnalyzeRlist$isoformFeatures$isoform_id[match(
                tolower(myCPCresults$id),
                tolower(switchAnalyzeRlist$isoformFeatures$isoform_id)
            )]

        ### subset
        myCPCresults <-
            myCPCresults[which(
                myCPCresults$id %in%
                    switchAnalyzeRlist$isoformFeatures$isoform_id
            ), ]
    }

    ### Add analysis to switchAnalyzeRlist
    if (TRUE) {
        matchVector <-
            match(switchAnalyzeRlist$isoformFeatures$isoform_id,
                  myCPCresults$id)

        switchAnalyzeRlist$isoformFeatures$codingPotentialValue <-
            myCPCresults$coding_probability[matchVector]
        switchAnalyzeRlist$isoformFeatures$codingPotential      <-
            myCPCresults$coding_probability[matchVector] >= codingCutoff


        n <- length(unique(myCPCresults$id))
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
                ! switchAnalyzeRlist$isoformFeatures$codingPotential
            )])

        ### Replace noncoding isoforms
        toModify <- which( switchAnalyzeRlist$orfAnalysis$isoform_id %in% nonCodingIsoforms)
        switchAnalyzeRlist$orfAnalysis[
            toModify,
            which( ! colnames(switchAnalyzeRlist$orfAnalysis) %in% c('isoform_id','orf_origin'))
        ] <- NA
        if( ! is.null(switchAnalyzeRlist$orfAnalysis$orf_origin)) {
            switchAnalyzeRlist$orfAnalysis$orf_origin[toModify] <- 'Predicted'
        }


        ### Overwrite PTC
        switchAnalyzeRlist$isoformFeatures$PTC[which(
            switchAnalyzeRlist$isoformFeatures$isoform_id %in%
                nonCodingIsoforms
        )] <- NA
    }

    return(switchAnalyzeRlist)
}

analyzePFAM <- function(
    ### Core arguments
    switchAnalyzeRlist,
    pathToPFAMresultFile,
    
    ### Advanced arguments
    ignoreAfterBar = TRUE,
    ignoreAfterSpace = TRUE,
    ignoreAfterPeriod = FALSE,
    showProgress = TRUE,
    quiet = FALSE
) {
    ### Test input
    if(TRUE) {
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(
                'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
            )
        }
        if (is.null(switchAnalyzeRlist$orfAnalysis)) {
            stop('ORF needs to be analyzed. Please run \'addORFfromGTF()\' (and if nessesary \'analyzeNovelIsoformORF()\') and try again.')
        }

        # file
        if (class(pathToPFAMresultFile) != 'character') {
            stop(
                'The \'pathToPFAMresultFile\' argument must be a string pointing to the PFAM result file(s)'
            )
        }
        if ( ! all(sapply(pathToPFAMresultFile, file.exists)) ) {
            stop('The file(s) \'pathToPFAMresultFile\' points to does not exist')
        }
    }

    if (showProgress & !quiet) {
        progressBar <- 'text'
    } else {
        progressBar <- 'none'
    }

    ### Import result data
    if(TRUE) {
        ### Figure out file type
        if(TRUE) {
            fileVector <- readLines(pathToPFAMresultFile[1], n = 3)

            if (! length(fileVector) ) {
                stop('The file pointed to by \'pathToPFAMresultFile\' is empty')
            }

            tapsAnnotated <- any(stringr::str_detect(fileVector, '\\t'))
        }
        # tapsAnnotated

        ### Figure out header and lines to skip
        if(TRUE) {
            ### Test whether headers are included
            temp <-
                read.table(
                    file = pathToPFAMresultFile[1],
                    stringsAsFactors = FALSE,
                    fill = TRUE,
                    header = FALSE,
                    nrows = 1
                )
            headerIncluded <- grepl('^<seq|^seq', temp[1, 1])

            ### Figure out number of lines to skip
            if(   headerIncluded) {
                tmp2 <-
                    readLines(
                        con = pathToPFAMresultFile[1],
                        n = 50
                    )
                skipLine <- which(grepl('^<seq|^seq', tmp2))[1]

                temp3 <-
                    read.table(
                        file = pathToPFAMresultFile[1],
                        stringsAsFactors = FALSE,
                        fill = TRUE,
                        header = FALSE,
                        nrows = 1,
                        skip = skipLine,
                    )
            }
            if( ! headerIncluded) {
                skipLine <- 0
            }
        }
        # skipLine

        ### Read in file(s)
        myPfamResult <- unique(
            do.call(
                rbind,
                plyr::llply(
                    pathToPFAMresultFile,
                    .fun = function(
                        aFile
                    ) {
                        # aFile <- pathToPFAMresultFile[1]

                        ### Read in file using the filte type info extracted above
                        if(TRUE) {
                            if( tapsAnnotated ) {
                                localPfamRes <- read.table(
                                    file = aFile,
                                    stringsAsFactors = FALSE,
                                    fill = TRUE,
                                    header = FALSE,
                                    col.names = 1:16,
                                    sep = '\t',
                                    skip = skipLine
                                )
                            } else {
                                localPfamRes <- read.table(
                                    file = aFile,
                                    stringsAsFactors = FALSE,
                                    fill = TRUE,
                                    header = FALSE,
                                    col.names = 1:16,
                                    sep = '', # one or more white spaces
                                    skip = skipLine
                                )
                            }

                        }

                        ### Determine if active sites are predicted
                        activeResidueIncluded <- any(
                            stringr::str_detect(
                                string = localPfamRes[,ncol(localPfamRes)],
                                pattern = 'predicted_active_site'
                            ),
                            na.rm = TRUE
                        )

                        ### Test shifted values
                        testShiftedValues <- function(aDF) {
                            try1 <- unique(c(
                                which(is.na( aDF[,14]     )),              # via "significant" column (will either be NA or a clan indication)
                                which(       aDF[,14] != 1 ),              # via "significant" column (will either be NA or a clan indication)
                                which( ! stringr::str_detect(aDF[,6], '^PF|^PB') )  # via pfam_hmm id which should start with PF or PB
                            ))

                            return(try1)
                        }
                        try1problems <- testShiftedValues(localPfamRes)

                        ### If any NAs - try reading it with read_fwf
                        if( length(try1problems) ) {
                            ### Fix the mistake in the fixed width file
                            suppressWarnings(
                                fileVector <- readLines(aFile, encoding = 'UTF-8')
                            )
                            fileVector <- gsub('Coiled-coil',' Coiled', fileVector)

                            ### FOr read_fwf to recognice it as raw text it must have new lines included
                            fileVector <- stringr::str_c(fileVector, '\n')

                            suppressMessages(
                                suppressWarnings(
                                    localPfamRes2 <- readr::read_fwf(
                                        file = fileVector,
                                        col_positions = readr::fwf_empty(
                                            file = fileVector,
                                            col_names = paste0(
                                                'X',
                                                1:(15+ activeResidueIncluded)
                                            ),
                                            skip = skipLine,
                                            comment = '#',
                                            n = 1000L
                                        ),
                                        skip = skipLine,
                                        comment = '#'
                                    )
                                )
                            )
                            localPfamRes2 <- as.data.frame(localPfamRes2)

                            ### Test shifted values via "significant" column (will either be NA or a clan indication)
                            try2problems <- testShiftedValues(localPfamRes2)

                            ### overwrite if fewer problems
                            if( length(try2problems) < length(try1problems) ) {
                                localPfamRes <- localPfamRes2
                            }
                        }

                        return(localPfamRes)
                    }
                )
            )
        )

        ### read in pfam resut result
        if (nrow(myPfamResult) == 0) {
            stop('The file pointed to by \'pathToPFAMresultFile\' is empty')
        }

        ### Identify they style of PFAM file(s) is supplied
        if(TRUE) {
            ### Define Colnames for pfam file types
            if(TRUE) {
                oldColnames <-
                    c(
                        'seq_id',
                        'alignment_start',
                        'alignment_end',
                        'envelope_start',
                        'envelope_end',
                        'hmm_acc',
                        'hmm_name',
                        'type', # 8
                        'hmm_start',
                        'hmm_end',
                        'hmm_length',
                        'bit_score',
                        'E_value',
                        'significant', # 14
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
                        'hmm_start', # 8
                        'hmm_end',
                        'hmm_length',
                        'bit_score',
                        'Individual_E_value',
                        'Conditional_E_value',
                        'significant', # 14
                        'outcompeted',
                        'clan'
                    )
            }

            ### Define style via clan column
            oldStyle <- any(grepl('^CL', myPfamResult[,15]))
            if(ncol(myPfamResult) == 16 ) {
                newStyle <- any(grepl('^CL', myPfamResult[,16]))
            } else {
                newStyle <- FALSE
            }

            ### Old style
            if ( oldStyle & !newStyle ) {
                colnames(myPfamResult) <- oldColnames[1:ncol(myPfamResult)]

                if( 'residue' %in% colnames(myPfamResult) ) {
                    myPfamResult$residue[which(myPfamResult$residue == '')] <-
                        NA
                }
            }

            ### New style
            if ( ! oldStyle & newStyle ) {
                colnames(myPfamResult) <- newColNames[1:ncol(myPfamResult)]
                myPfamResult$clan[which(myPfamResult$clan == '')] <- NA

            }

            ### Other combinations
            if( oldStyle == newStyle ) {
                stop('The file(s) supplied is not recogniced as a pfam output.')
            }
        }

        ### Revert the coiled name
        if('type' %in% colnames(myPfamResult)) {
            myPfamResult$type <- gsub('^Coiled$','Coiled_coil', myPfamResult$type)
        }

    }

    ### Sanity check that it is a PFAM result file
    if (TRUE) {
        ### Test expected clumns
        test1 <-
            ncol(myPfamResult) == 15 |
            ncol(myPfamResult) == 16 # the output have 15 or 16 collumns depending on whther active sites are predicted

        if( any(is.na(myPfamResult$hmm_acc)) ) {
            naIndex <- which(
                is.na(myPfamResult$hmm_acc) | is.na(myPfamResult$type)
            )

            problemSize <- length(naIndex)

            exampleProblems <- sample(
                myPfamResult$hmm_name[naIndex],
                size = min(3, problemSize)
            )

            warning(
                paste(
                    'There are',
                    problemSize,
                    'enteries in the pfam results files where there are missing values.',
                    '\nMost likely these will not affect the analysis done with IsoformSwitchAnalyzeR',
                    '\nbut we recomend the user to figure out why there are missing data in the first place.',
                    '\nExamples of hmm_names with missing data are:',
                    paste(exampleProblems, collapse = ', '),
                    sep = ' '
                )
            )

            nonNaIndex <- which(
                ! is.na(myPfamResult$hmm_acc)
            )
        } else {
            nonNaIndex <- 1:nrow(myPfamResult)
        }
        test2 <-
            all(grepl(
                pattern = '^PF|^PB' ,
                myPfamResult$hmm_acc[nonNaIndex],
                ignore.case = FALSE
            ))                # All pfam hmm starts with PF or PB

        if (!all(test1, test2)) {
            stop('The file pointed to by \'pathToPFAMresultFile\' is not a PFAM result file')
        }

        ### Potentially remove '>'
        t1 <-   all( grepl( '^>', myPfamResult$seq_id ))
        t2 <- ! all( grepl( '^>', switchAnalyzeRlist$isoformFeatures$isoform_id ))
        if( t1 & t2 ) {
            myPfamResult$seq_id <- gsub('^>', '', myPfamResult$seq_id)
            warning('Removed the prefix ">" from all Pfam results since we suspect they are not supposed to be there.')
        }
        
        ### Fix names
        if( ignoreAfterBar | ignoreAfterSpace | ignoreAfterPeriod){
          myPfamResult$seq_id <- fixNames(
            nameVec = myPfamResult$seq_id,
            ignoreAfterBar = ignoreAfterBar,
            ignoreAfterSpace = ignoreAfterSpace,
            ignoreAfterPeriod = ignoreAfterPeriod
          )
        }

        ### test names
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

    ### Add Domain isotypes
    if(TRUE) {
        tmp <-
            myPfamResult %>%
            pfamAnalyzeR::augment_pfam() %>%
            pfamAnalyzeR::analyse_pfam_isotypes()

        myPfamResult$domain_isotype        <- tmp$domain_isotype
        myPfamResult$domain_isotype_simple <- tmp$domain_isotype_simple

    }

    ### Convert from AA coordinates to transcript and gnomic coordinats
    if (TRUE) {
        if (!quiet) {
            message('Converting AA coordinats to transcript and genomic coordinats...')
        }
        ### Remove unwanted columns
        myPfamResult$alignment_start <- NULL # KVS update 01/2023
        myPfamResult$alignment_end   <- NULL # KVS update 01/2023
        myPfamResult$hmm_start <- NULL
        myPfamResult$hmm_end <- NULL
        myPfamResult$hmm_length <- NULL

        colnames(myPfamResult)[which(
            grepl('envelope_', colnames(myPfamResult)) # KVS update 01/2023
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

        #myPfamResultDf$pfamStarExon <- NULL
        #myPfamResultDf$pfamEndExon <- NULL

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

    ### Repport numbers
    if(TRUE) {
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
    }

    return(switchAnalyzeRlist)
}

analyzeSignalP <- function(
    ### Core arguments
    switchAnalyzeRlist,
    pathToSignalPresultFile,
    
    ### Advanced arguments
    minSignalPeptideProbability = 0.5,
    ignoreAfterBar = TRUE,
    ignoreAfterSpace = TRUE,
    ignoreAfterPeriod = FALSE,
    quiet = FALSE
) {
    ### Test input
    if(TRUE) {
        if (is.null(switchAnalyzeRlist$orfAnalysis)) {
            stop('ORF needs to be analyzed. Please run \'addORFfromGTF()\' (and if nessesary \'analyzeNovelIsoformORF()\') and try again.')
        }

        # file
        if (class(pathToSignalPresultFile) != 'character') {
            stop(
                'The \'pathToSignalPresultFile\' argument must be a string pointing to the SignalP result file(s)'
            )
        }
        if ( ! all(sapply(pathToSignalPresultFile, file.exists)) ) {
            stop('The file(s) \'pathToSignalPresultFile\' points to does not exist')
        }

    }

    ### Obtain signalP result
    if (TRUE) {
        ### Figure out which version of signal-p was used
        fileType <- read.table(
            pathToSignalPresultFile[1],
            nrows = 1,
            fill = TRUE,
            header = FALSE,
            comment.char = '',
            sep='\t',
            stringsAsFactors = FALSE
        )

        ### Annotate version
        if( grepl('SignalP-5', fileType$V1[1]) ) {
            signalPversion = 'v5'
        } else if( grepl('SignalP-6', fileType$V1[1]) ) {
            signalPversion = 'v6'
        } else {
            signalPversion = 'v4'
        }

        ### SignalP6
        if( signalPversion == 'v6' ) {
            if( !any( sapply(unlist(fileType), function(x) {grepl('Eukarya|euk',x)} )) ){
                warning('It seems SignalP was run as Non-Eukaryote - was that on purpouse?')
            }

            signalP6colNames <- c(
                'isoform_id',
                'Prediction',
                'OTHER',
                'SP_Sec_SPI',
                'CS_Position'
            )

            ### Read in SignalP predictions
            if(TRUE) {
                singalPresults <- do.call(rbind, plyr::llply(
                    pathToSignalPresultFile,
                    .fun = function(
                        aFile
                    ) {
                        suppressWarnings(
                            readr::read_tsv(
                                file = aFile,
                                col_names = signalP6colNames,
                                col_types = readr::cols(
                                    isoform_id  = readr::col_character(),
                                    Prediction  = readr::col_character(),
                                    SP_Sec_SPI  = readr::col_double(),
                                    OTHER       = readr::col_double(),
                                    CS_Position = readr::col_character()
                                ),
                                comment = '#'
                            )
                        )
                        #read.table(
                        #    aFile,
                        #    header = FALSE,
                        #    stringsAsFactors = FALSE,
                        #    fill = TRUE,
                        #    col.names = signalP5colNames,
                        #    sep='\t'
                        #)
                    }
                ))

                colnames(singalPresults) <- gsub('\\.$','', colnames(singalPresults))
                colnames(singalPresults) <- gsub('\\.','_', colnames(singalPresults))

                singalPresults <- unique(singalPresults)
            }

            ### Sanity check that it is a SignalIP result file
            if(TRUE) {
                if(nrow(singalPresults) == 0) {
                    stop('The result file(s) seems to be empty')
                }

                t1 <- ! is.character(singalPresults$isoform_id)
                t2 <- ! is.character(singalPresults$Prediction)
                t3 <- ! is.numeric(singalPresults$SP_Sec_SPI)
                t4 <- ! is.numeric(singalPresults$OTHER)
                t5 <- ! any(c('OTHER','NO_SP','SP') %in% singalPresults$Prediction) # add "Other"

                if( any( c(t1,t2,t3,t4,t5))) {
                    stop('The pathToSignalPresultFile does not seam to be the result of a SignalP 6 analysis')
                }
                
                if( ignoreAfterBar | ignoreAfterSpace | ignoreAfterPeriod){
                  singalPresults$transcript_id <- fixNames(
                    nameVec = singalPresults$transcript_id,
                    ignoreAfterBar = ignoreAfterBar,
                    ignoreAfterSpace = ignoreAfterSpace,
                    ignoreAfterPeriod = ignoreAfterPeriod
                  )
                }

                if( ! any( singalPresults$isoform_id %in% switchAnalyzeRlist$isoformFeatures$isoform_id) ) {
                    stop('The pathToSignalPresultFile does not contain result of isoforms analyzed in the switchAnalyzeRlist')
                }
            }

            ### Reduce to features of interest
            if(TRUE) {
                ### With signal
                singalPresults <- singalPresults[which(
                    singalPresults$Prediction == 'SP'
                ),]

                ### Also subset on those with cleaveage site annotated
                singalPresults <- singalPresults[which(
                    grepl('CS pos:', singalPresults$CS_Position)
                ),]

                ### Remove fragment
                singalPresults <- singalPresults[which(
                    ! grepl('CS pos: \\?. Probable protein fragment', singalPresults$CS_Position)
                ),]

                ### Test prob
                toLargeProb <- unique( c(
                    which( singalPresults$SP_Sec_SPI > 1),
                    which( singalPresults$OTHER > 1)
                ))
                if(length(toLargeProb)) {
                    warning(
                        paste(
                            'Please note that', length(toLargeProb),
                            'SignalP results had probabilities ',
                            '\nlarger than 1 - which probabilities cannot be.',
                            '\nThese enteries were removed.',
                            sep = ' '
                        )
                    )
                    singalPresults <- singalPresults[setdiff( 1:nrow(singalPresults), toLargeProb),]
                }


                # Analyzed with ORF
                singalPresults <- singalPresults[which(
                    singalPresults$isoform_id %in% switchAnalyzeRlist$orfAnalysis$isoform_id[which(
                        !is.na(switchAnalyzeRlist$orfAnalysis$orfTransciptStart)
                    )]
                ),]

                if( nrow(singalPresults) == 0) {
                    warning('There were no signal peptides in the SignalP result file. Returning switchAnalyzeRlist without this analysis.')
                    return(switchAnalyzeRlist)
                }
            }

            ### Massage
            if(TRUE) {
                singalPresults$aa_removed <- sapply(
                    strsplit(singalPresults$CS_Position, split = '\\.|-'),
                    function(x) {
                        as.integer( gsub('^CS pos: ', '', x[1]) )
                    }
                )

                singalPresults$Prediction <- NULL
            }
        }

        ### SignalP5
        if( signalPversion == 'v5' ) {
            if( !any( sapply(unlist(fileType), function(x) {grepl('Eukarya|euk',x)} )) ){
                warning('It seems SignalP was run as Non-Eukaryote - was that on purpouse?')
            }

            signalP5colNames <- c(
                'isoform_id',
                'Prediction',
                'SP_Sec_SPI',
                'OTHER',
                'CS_Position'
            )

            ### Read in SignalP predictions
            if(TRUE) {
                singalPresults <- do.call(rbind, plyr::llply(
                    pathToSignalPresultFile,
                    .fun = function(
                        aFile
                    ) {
                        suppressWarnings(
                            readr::read_tsv(
                                file = aFile,
                                col_names = signalP5colNames,
                                col_types = readr::cols(
                                    isoform_id  = readr::col_character(),
                                    Prediction  = readr::col_character(),
                                    SP_Sec_SPI  = readr::col_double(),
                                    OTHER       = readr::col_double(),
                                    CS_Position = readr::col_character()
                                ),
                                comment = '#'
                            )
                        )
                        #read.table(
                        #    aFile,
                        #    header = FALSE,
                        #    stringsAsFactors = FALSE,
                        #    fill = TRUE,
                        #    col.names = signalP5colNames,
                        #    sep='\t'
                        #)
                    }
                ))

                colnames(singalPresults) <- gsub('\\.$','', colnames(singalPresults))
                colnames(singalPresults) <- gsub('\\.','_', colnames(singalPresults))

                singalPresults <- unique(singalPresults)
            }

            ### Sanity check that it is a SignalIP result file
            if(TRUE) {
                if(nrow(singalPresults) == 0) {
                    stop('The result file(s) seems to be empty')
                }

                t1 <- ! is.character(singalPresults$isoform_id)
                t2 <- ! is.character(singalPresults$Prediction)
                t3 <- ! is.numeric(singalPresults$SP_Sec_SPI)
                t4 <- ! is.numeric(singalPresults$OTHER)
                t5 <- ! all(singalPresults$Prediction %in% c('OTHER','SP(Sec/SPI)','LIPO(Sec/SPII)','TAT(Tat/SPI)'))

                if( any( c(t1,t2,t3,t4,t5))) {
                    stop('The pathToSignalPresultFile does not seam to be the result of a SignalP 5 analysis')
                }
                
                if( ignoreAfterBar | ignoreAfterSpace | ignoreAfterPeriod){
                  singalPresults$transcript_id <- fixNames(
                    nameVec = singalPresults$transcript_id,
                    ignoreAfterBar = ignoreAfterBar,
                    ignoreAfterSpace = ignoreAfterSpace,
                    ignoreAfterPeriod = ignoreAfterPeriod
                  )
                }

                if( ! any( singalPresults$isoform_id %in% switchAnalyzeRlist$isoformFeatures$isoform_id) ) {
                    stop('The pathToSignalPresultFile does not contain result of isoforms analyzed in the switchAnalyzeRlist')
                }
            }

            ### Reduce to features of interest
            if(TRUE) {
                peptideCols <- c('SP_Sec_SPI','TAT_Tat_SPI','LIPO_Sec_SPII')

                ### With signal
                singalPresults <- singalPresults[which(
                    apply(
                        X = singalPresults[,na.omit(match(peptideCols, colnames(singalPresults))), drop=FALSE],
                        MARGIN = 1,
                        FUN = function(x) { any(x >= minSignalPeptideProbability)}

                    )
                ),]

                ### Also subset on those with cleaveage site annotated
                singalPresults <- singalPresults[which(
                    grepl('CS pos:', singalPresults$CS_Position)
                ),]

                ### Remove fragment
                singalPresults <- singalPresults[which(
                    ! grepl('CS pos: \\?. Probable protein fragment', singalPresults$CS_Position)
                ),]

                ### Test prob
                toLargeProb <- unique( c(
                    which( singalPresults[,na.omit(match(peptideCols, colnames(singalPresults))), drop=TRUE] > 1),
                    which( singalPresults$OTHER > 1)
                ))
                if(length(toLargeProb)) {
                    warning(
                        paste(
                            'Please note that', length(toLargeProb),
                            'SignalP results had probabilities ',
                            '\nlarger than 1 - which probabilities cannot be.',
                            '\nThese enteries were removed.',
                            sep = ' '
                        )
                    )
                    singalPresults <- singalPresults[setdiff( 1:nrow(singalPresults), toLargeProb),]
                }


                # Analyzed with ORF
                singalPresults <- singalPresults[which(
                    singalPresults$isoform_id %in% switchAnalyzeRlist$orfAnalysis$isoform_id[which(
                        !is.na(switchAnalyzeRlist$orfAnalysis$orfTransciptStart)
                    )]
                ),]

                if( nrow(singalPresults) == 0) {
                    warning('There were no signal peptides in the SignalP result file. Returning switchAnalyzeRlist without this analysis.')
                    return(switchAnalyzeRlist)
                }
            }

            ### Massage
            if(TRUE) {
                singalPresults$aa_removed <- sapply(
                    strsplit(singalPresults$CS_Position, split = '\\.|-'),
                    function(x) {
                        as.integer( gsub('^CS pos: ', '', x[1]) )
                    }
                )

                singalPresults$Prediction <- NULL
            }
        }

        ### SignalP4
        if( signalPversion == 'v4' ) {
            ### Old file
            singalPresults <- do.call(rbind, plyr::llply(
                pathToSignalPresultFile,
                .fun = function(
                    aFile
                ) {
                    read.table(
                        aFile,
                        header = FALSE,
                        stringsAsFactors = FALSE,
                        fill = TRUE,
                        col.names = paste('V', 1:13, sep = '')
                    )
                }
            ))
            if(nrow(singalPresults) == 0) {
                stop('The result file(s) seems to be empty')
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
                    stop('No signial peptides were found in the isoforms analyzed')
                }

                singalPresults <- unique(singalPresults)


                t1 <- all(singalPresults$V3 == 'Cleavage')
                t2 <- all(singalPresults$V4 == 'site')
                t3 <- all(singalPresults$V5 == 'between')
                t4 <- all(singalPresults$V6 == 'pos.')
                t5 <- all(grepl('^Name=', singalPresults$V1))
                t6 <- all(grepl('^Networks=', singalPresults$V13))

                if (!all(t1, t2, t3, t4, t5, t6)) {
                    stop(
                        'The file(s) pointed to by \'pathToSignalPresultFile\' is not a SignalIP summary file'
                    )
                }
            }

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
            if( ignoreAfterBar | ignoreAfterSpace | ignoreAfterPeriod){
              singalPresults$transcript_id <- fixNames(
                nameVec = singalPresults$transcript_id,
                ignoreAfterBar = ignoreAfterBar,
                ignoreAfterSpace = ignoreAfterSpace,
                ignoreAfterPeriod = ignoreAfterPeriod
              )
            }

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

            colnames(singalPresults)[1] <- 'isoform_id'

        }

    }

    ### Convert from AA coordinats to transcript and genomic coordinats
    if (TRUE) {
        ### Convert cleavage site to transcript coordinats
        # df with orf start
        orfStartDF <-
            unique(as.data.frame(switchAnalyzeRlist$orfAnalysis[,c(
                'isoform_id', 'orfTransciptStart'
            )]))
        # calculate start position
        singalPresults$transcriptClevageAfter <-
            (singalPresults$aa_removed  * 3) + orfStartDF[match(
                x = singalPresults$isoform_id,
                table = orfStartDF$isoform_id
            ), 2] - 1

        ### convert from transcript to genomic coordinats
        # extract exon data
        myExons <-
            as.data.frame(switchAnalyzeRlist$exons[which(
                switchAnalyzeRlist$exons$isoform_id %in%
                    singalPresults$isoform_id
            ),])
        myExonsSplit <- split(myExons, f = myExons$isoform_id)

        # calculate genomic coordinat
        singalPresults <-
            plyr::ddply(
                singalPresults,
                .variables = 'isoform_id',
                .fun = function(aDF) {
                    localExons <- myExonsSplit[[aDF$isoform_id[1]]]

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

analyzeNetSurfP3 <- function(
    ### Advanced arguments
    switchAnalyzeRlist,
    pathToNetSurfP3resultFile,
    
    ### Advanced arguments
    smoothingWindowSize = 5,
    probabilityCutoff = 0.5,
    minIdrSize = 30,
    ignoreAfterBar = TRUE,
    ignoreAfterSpace = TRUE,
    ignoreAfterPeriod = FALSE,
    showProgress = TRUE,
    quiet = FALSE
) {
  ### Test input
  if(TRUE) {
    if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
      stop(
        'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
      )
    }
    if (is.null(switchAnalyzeRlist$orfAnalysis)) {
      stop('ORF needs to be analyzed. Please run \'addORFfromGTF()\' (and if nessesary \'analyzeNovelIsoformORF()\') and try again.')
    }
    
    # file
    if (class(pathToNetSurfP3resultFile) != 'character') {
      stop(
        'The \'pathToNetSurfP3resultFile\' argument must be a string pointing to the NetSurfP3 result file'
      )
    }
    if (! all( file.exists(pathToNetSurfP3resultFile)) ) {
      stop('(At least one of) the file(s) \'pathToNetSurfP3resultFile\' points to does not exist')
    }
    
    if( smoothingWindowSize %% 2 != 1 | !is(smoothingWindowSize, 'numeric') ) {
      stop('The \'smoothingWindowSize\' argument must be an odd integer')
    }
  }
  
  if (showProgress & !quiet) {
    progressBar <- 'text'
  } else {
    progressBar <- 'none'
  }
  
  ### Read result file
  if(TRUE) {
    if (!quiet) {
      message('Step 1 of 3: Reading results into R...')
    }
    
    ### Read in file
    suppressWarnings(
      netSurf <- do.call(rbind, plyr::llply(
        pathToNetSurfP3resultFile,
        .fun = function(
    aFile
        ) {
          read_csv(
            file = aFile,
            col_names = TRUE,
            col_types = cols_only(
              id = col_character(),
              n = col_integer(),
              disorder = col_double()
            ),
            progress = showProgress & !quiet
          ) %>%
            dplyr::mutate(id = sub("^>", "", id))
        }
      ))
    
    )
    ### Fix names
    if( ignoreAfterBar | ignoreAfterSpace | ignoreAfterPeriod){
      netSurf$id <- fixNames(
        nameVec = netSurf$id,
        ignoreAfterBar = ignoreAfterBar,
        ignoreAfterSpace = ignoreAfterSpace,
        ignoreAfterPeriod = ignoreAfterPeriod
      )
    }
    
    ### Sanity check
    if( ! any(netSurf$id %in% switchAnalyzeRlist$isoformFeatures$isoform_id )) {
      stop('The \'pathToNetSurfP3resultFile\' does not appear to contain results for the isoforms stored in the switchAnalyzeRlist...')
    }
    
    ### Subset to relevant features
    netSurf <- netSurf[which(
      netSurf$id %in%
        switchAnalyzeRlist$orfAnalysis$isoform_id[which(
          !is.na(switchAnalyzeRlist$orfAnalysis$orfTransciptStart)
        )]
    ),]
    
    netSurf <- unique(netSurf)
  }
  
  ### Reduce to those with IDR
  if(TRUE) {
    if (!quiet) {
      message('Step 2 of 3: Analyzing data to extract IDRs...')
    }
    netSurf$idDis <- as.integer( netSurf$disorder > probabilityCutoff )
    
    
    disRle <- RleList(
      split(round(netSurf$disorder, digits = 3), netSurf$id)
    )
    
    ### Apply spliding window
    disRle <- runmean(disRle, k=smoothingWindowSize, endrule = 'drop')
    nRemovedByDrop <- (smoothingWindowSize-1) / 2
    
    ### Convert to binary
    disRle <- disRle > probabilityCutoff
    
    ### Loop over and extract result
    disRes <- plyr::ldply(disRle, .progress = progressBar, function(localRle) {
      ### Extract start and stop
      rleDf <- data.frame(
        classification=localRle@values,
        length=localRle@lengths,
        orf_aa_end=cumsum(localRle@lengths) + nRemovedByDrop
      )
      rleDf$orf_aa_start <- rleDf$orf_aa_end - rleDf$length + 1 + nRemovedByDrop
      
      ### Subset to disordered of length X
      rleDf <- rleDf[which(
        rleDf$classification &
          rleDf$length >= minIdrSize
      ),]
      
      rleDf$classification <- NULL
      
      return(rleDf)
    })
    colnames(disRes)[1] <- 'isoform_id'
    
    disRes <- disRes[,c('isoform_id','orf_aa_start','orf_aa_end','length')]
    
    ### Add type
    disRes$idr_type <- 'IDR'
  }
  
  ### Convert from AA coordinates to transcript and genomic coordinates
  if (TRUE) {
    if (!quiet) {
      message('Step 3 of 3: Converting AA coordinates to transcript and genomic coordinates...')
    }
    
    ### convert from codons to transcript position
    orfStartDF <-
      unique(
        switchAnalyzeRlist$orfAnalysis[
          which( !is.na(switchAnalyzeRlist$orfAnalysis$orfTransciptStart)),
          c('isoform_id', 'orfTransciptStart')
        ]
      )
    disRes$transcriptStart <-
      (disRes$orf_aa_start  * 3 - 2) +
      orfStartDF[
        match(
          x = disRes$isoform_id,
          table = orfStartDF$isoform_id
        ),
        2] - 1
    disRes$transcriptEnd <-
      (disRes$orf_aa_end * 3) +
      orfStartDF[
        match(
          x = disRes$isoform_id,
          table = orfStartDF$isoform_id
        ),
        2] - 1
    
    ### convert from transcript to genomic coordinates
    # extract exon data
    myExons <-
      as.data.frame(switchAnalyzeRlist$exons[which(
        switchAnalyzeRlist$exons$isoform_id %in% disRes$isoform_id
      ), ])
    myExonsSplit <- split(myExons, f = myExons$isoform_id)
    
    # loop over the individual transcripts and extract the genomic coordiants of the domain and also for the active residues (takes 2 min for 17000 rows)
    disResDf <-
      plyr::ddply(
        disRes,
        .progress = progressBar,
        .variables = 'isoform_id',
        .fun = function(aDF) {
          # aDF <- disRes[which(disRes$isoform_id == 'uc001isa.1'),]
          
          transcriptId <- aDF$isoform_id[1]
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
            
            orfPosList[[as.character(j)]] <- domainInfo
            
          }
          orfPosDf <- do.call(rbind, orfPosList)
          
          return(cbind(aDF, orfPosDf))
        }
      )
    
    colnames(disResDf) <- gsub(
      'pfam',
      'idr',
      colnames(disResDf)
    )
    
  }
  
  ### Add analysis to switchAnalyzeRlist
  if (TRUE) {
    # sort
    disResDf <-
      disResDf[order(
        disResDf$isoform_id,
        disResDf$transcriptStart
      ), ]
    
    #disResDf$idrStarExon <- NULL
    #disResDf$idrEndExon <- NULL
    
    # add the results to the switchAnalyzeRlist object
    switchAnalyzeRlist$idrAnalysis <- disResDf
    
    # add indication to transcriptDf
    switchAnalyzeRlist$isoformFeatures$idr_identified <- 'no'
    switchAnalyzeRlist$isoformFeatures$idr_identified[which(
      is.na(switchAnalyzeRlist$isoformFeatures$PTC)
    )] <- NA # sets NA for those not analyzed
    
    switchAnalyzeRlist$isoformFeatures$idr_identified[which(
      switchAnalyzeRlist$isoformFeatures$isoform_id %in%
        disResDf$isoform_id
    )] <- 'yes'
  }
  
  n <- length(unique(disResDf$isoform_id))
  p <-
    round(n / length(unique(
      switchAnalyzeRlist$isoformFeatures$isoform_id
    )) * 100, digits = 2)
  
  if (!quiet) {
    message(paste(
      'Added IDR information to ',
      n,
      ' (',
      p,
      '%) transcripts',
      sep = ''
    ))
  }
  return(switchAnalyzeRlist)
}

analyzeIUPred2A <- function(
    ### Core arguments
    switchAnalyzeRlist,
    pathToIUPred2AresultFile,
    
    ### Advanced arguments
    smoothingWindowSize = 11,
    probabilityCutoff = 0.5,
    minIdrSize = 30,
    annotateBindingSites = TRUE,
    minIdrBindingSize = 15,
    minIdrBindingOverlapFrac = 0.8,
    ignoreAfterBar = TRUE,
    ignoreAfterSpace = TRUE,
    ignoreAfterPeriod = FALSE,
    showProgress = TRUE,
    quiet = FALSE
) {
    ### Test input
    if(TRUE) {
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(
                'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
            )
        }
        if (is.null(switchAnalyzeRlist$orfAnalysis)) {
            stop('ORF needs to be analyzed. Please run \'addORFfromGTF()\' (and if nessesary \'analyzeNovelIsoformORF()\') and try again.')
        }

        # file
        if (class(pathToIUPred2AresultFile) != 'character') {
            stop(
                'The \'pathToIUPred2AresultFile\' argument must be a string pointing to the IUPred2A result file'
            )
        }
        if (! all( file.exists(pathToIUPred2AresultFile)) ) {
            stop('(At least on of) the file(s) \'pathToIUPred2AresultFile\' points to does not exist')
        }

        if( smoothingWindowSize %% 2 != 1 | !is(smoothingWindowSize, 'numeric') ) {
            stop('The \'smoothingWindowSize\' argument must be an odd integer')
        }

        totalAnalysisNr <- 3 + as.integer(annotateBindingSites)
        analysisNr <- 1
    }

    if (showProgress & !quiet) {
        progressBar <- 'text'
    } else {
        progressBar <- 'none'
    }

    ### Read result file
    if(TRUE) {
        if (!quiet) {
            message(
                paste(
                    'Step', analysisNr, 'of', totalAnalysisNr,
                    ': Reading results into R...'
                )
            )
        }

        ### Read in file
        suppressWarnings(
            iupred2a <- do.call(rbind, plyr::llply(
                pathToIUPred2AresultFile,
                .fun = function(
                    aFile
                ) {
                    read.table(
                        file = aFile,
                        header = FALSE,
                        sep = '\t',
                        stringsAsFactors = FALSE,
                        strip.white = TRUE,
                        comment.char = '#',
                        blank.lines.skip = TRUE,
                        fill = TRUE
                    )
                }
            ))
        )

        ### Massage into tidy format
        if(TRUE) {
            ### Fix fails in reading in files
            if(TRUE) {
                ### Remove comments
                iupred2a <- iupred2a[which(
                    ! grepl('^#', iupred2a$V1)
                ),]
                ### Remove blank lines
                iupred2a <- iupred2a[which(
                    ! iupred2a$V1 == ''
                ),]
            }

            myNames <- gsub(
                '^>',
                '',
                iupred2a$V1[which(is.na(iupred2a$V3))]
            ) 
            
            ### Fix names
            if( ignoreAfterBar | ignoreAfterSpace | ignoreAfterPeriod){
                myNames <- fixNames(
                nameVec = myNames,
                ignoreAfterBar = ignoreAfterBar,
                ignoreAfterSpace = ignoreAfterSpace,
                ignoreAfterPeriod = ignoreAfterPeriod
              )
            }

            ### Sanity check
            if( ! any( myNames %in% switchAnalyzeRlist$isoformFeatures$isoform_id )) {
                stop('The \'pathToIUPred2AresultFile\' does not appear to contain results for the isoforms stored in the switchAnalyzeRlist...')
            }

            if( annotateBindingSites ) {
                if( ncol( iupred2a) != 4 ) {
                    stop(
                        paste(
                            'It appears the file pointed to by by pathToIUPred2AresultFile',
                            'does not contain the ANCHOR2 Predictions.',
                            'Please ensure the "context-dependent predictions" is turned on and set to ANCHOR2'
                        )
                    )
                }
            }

            ### Create id vector
            myIndex <- c(which(grepl('^>', iupred2a$V1)), nrow(iupred2a)+1)
            newIdDf <- data.frame(
                start = myIndex[-length(myIndex)],
                end = myIndex[-1] -1,
                newId = myNames,
                stringsAsFactors = FALSE
            )
            newIdDf$length <- newIdDf$end - newIdDf$start + 1

            iupred2a$newId <- rep(
                x = newIdDf$newId,
                times = newIdDf$length
            )
            if( any(is.na(iupred2a$newId)) ) {
                stop('Something went wrong in the convertion of to table format. Please contact developers with reproducible example')
            }

            ### Extract the parts of of interest
            iupred2a <- iupred2a[which( ! grepl('^>', iupred2a$V1)),]
            iupred2a <- iupred2a[,c(5,1:4)]

            colnames(iupred2a) <- c('isoform_id','position','residue','iupred_score','anchor_score')

            iupred2a$position <- as.integer(iupred2a$position)

            ### Sanity checks
            if(TRUE) {
                problemsFound <- FALSE
                if( ! is(iupred2a$iupred_score, 'numeric') ) {
                    problemsFound <- TRUE
                }
                if( ! is(iupred2a$anchor_score, 'numeric') ) {
                    problemsFound <- TRUE
                }
                if( ! all( iupred2a$residue %in% Biostrings::AA_ALPHABET ) ){
                    problemsFound <- TRUE
                }
                if(problemsFound) {
                    stop('The file(s) provided to "pathToIUPred2AresultFile" does not appear to be the result file from IUPred2A')
                }
            }
        }

        ### Subset to relecant features
        iupred2a <- iupred2a[which(
            iupred2a$isoform_id %in%
                switchAnalyzeRlist$orfAnalysis$isoform_id[which(
                    !is.na(switchAnalyzeRlist$orfAnalysis$orfTransciptStart)
                )]
        ),]

        iupred2a <- unique(iupred2a)

        analysisNr <- analysisNr + 1
    }

    ### Extract IDR regions of interest
    if(TRUE) {
        if (!quiet) {
            message(
                paste(
                    'Step', analysisNr, 'of', totalAnalysisNr,
                    ': Extracting regions of interest...'
                )
            )
        }

        ### Helper function
        if(TRUE) {
            extractRegionsAboveCutoff <- function(
                scoreList,
                scoreCutoff,
                minSize,
                smoothingWindowSize
            ) {
                scoreRle <- RleList( scoreList )

                ### Apply spliding window
                scoreRle <- runmean(scoreRle, k=smoothingWindowSize, endrule = 'drop')
                nRemovedByDrop <- (smoothingWindowSize-1) / 2

                ### Convert to binary
                scoreRle <- scoreRle > scoreCutoff

                ### Loop over and extract result
                regionRes <- plyr::ldply(scoreRle, .progress = progressBar, function(localRle) {
                    ### Extract start and stop
                    rleDf <- data.frame(
                        classification=localRle@values,
                        length=localRle@lengths,
                        orf_aa_end=cumsum(localRle@lengths) + nRemovedByDrop
                    )
                    rleDf$orf_aa_start <- rleDf$orf_aa_end - rleDf$length + 1 + nRemovedByDrop

                    ### Subset to disordered of length X
                    rleDf <- rleDf[which(
                        rleDf$classification &
                            rleDf$length >= minSize
                    ),]

                    rleDf$classification <- NULL

                    return(rleDf)
                })
                colnames(regionRes)[1] <- 'isoform_id'

                regionRes <- regionRes[,c('isoform_id','orf_aa_start','orf_aa_end','length')]
                return(regionRes)
            }
        }

        ### Extract IUPred2 IDR regions
        if(TRUE) {
            disRes <- extractRegionsAboveCutoff(
                scoreList = split(
                    iupred2a$iupred_score,
                    iupred2a$isoform_id
                ),
                scoreCutoff = probabilityCutoff,
                minSize = minIdrSize,
                smoothingWindowSize = smoothingWindowSize
            )

            disRes$idr_type <- 'IDR'
        }

        ### Extract ANCHOR2 IDR binding regions
        if( annotateBindingSites ) {
            bindRes <- extractRegionsAboveCutoff(
                scoreList = split(
                    iupred2a$anchor_score,
                    iupred2a$isoform_id
                ),
                scoreCutoff = probabilityCutoff,
                minSize = minIdrBindingSize,
                smoothingWindowSize = smoothingWindowSize
            )
        }

        analysisNr <- analysisNr + 1
    }

    ### Integrate the prediction of disordered regions and binding sites
    if( annotateBindingSites ) {
        if (!quiet) {
            message(
                paste(
                    'Step', analysisNr, 'of', totalAnalysisNr,
                    ': Integrating IDR with binding site predictions...'
                )
            )
        }

        disRes$rowId <- 1:nrow(disRes)

        ### Convert to GRanges
        disResGr <- GRanges(
            disRes$isoform_id,
            IRanges(
                disRes$orf_aa_start,
                disRes$orf_aa_end
            ),
            rowId = disRes$rowId
        )
        bindResGr <- GRanges(
            bindRes$isoform_id,
            IRanges(
                bindRes$orf_aa_start,
                bindRes$orf_aa_end
            )
        )

        ### Find overlap
        suppressWarnings(
            myOverlap <- findOverlaps(
                query   = disResGr,
                subject = bindResGr
            )
        )

        ### Intersect overlapping
        overlapIntersect <- suppressWarnings(
            pintersect(
                disResGr [queryHits(   myOverlap )],
                bindResGr[subjectHits( myOverlap )]
            )
        )
        ### Calculate fraction overlap (as fraction of binding region size)
        overlapIntersect$overlapFrac <- round(
            width(overlapIntersect) /
                width(bindResGr[subjectHits(myOverlap)]),
            digits = 3
        )

        ### Extract summary
        overlapFracList <- split(
            overlapIntersect$overlapFrac,
            overlapIntersect$rowId
        )
        overlapCount <- sapply(
            overlapFracList,
            length
        )
        overlapSize <- sapply(
            overlapFracList,
            max
        )

        ### Annotate overlap
        disRes$nr_idr_binding_region_overlapping <- NA
        overlapIndex <- match( names(overlapCount), disRes$rowId)
        disRes$nr_idr_binding_region_overlapping[ overlapIndex ] <- overlapCount

        disRes$max_fraction_of_idr_binding_region_overlapping <- NA
        overlapIndex <- match( names(overlapSize), disRes$rowId)
        disRes$max_fraction_of_idr_binding_region_overlapping[ overlapIndex ] <- overlapSize

        ### Chage type
        disRes$idr_type[which(
            disRes$max_fraction_of_idr_binding_region_overlapping > minIdrBindingOverlapFrac
        )] <- 'IDR_w_binding_region'

        ### Remove IDR id
        disRes$rowId <- NULL

        analysisNr <- analysisNr + 1
    }

    ### Convert from AA coordinats to transcript and genomic coordinats
    if (TRUE) {
        if (!quiet) {
            message(
                paste(
                    'Step', analysisNr, 'of', totalAnalysisNr,
                    ': Converting AA coordinats to transcript and genomic coordinats...'
                )
            )
        }

        ### convert from codons to transcript position
        orfStartDF <-
            unique(
                switchAnalyzeRlist$orfAnalysis[
                    which( !is.na(switchAnalyzeRlist$orfAnalysis$orfTransciptStart)),
                    c('isoform_id', 'orfTransciptStart')
                    ]
            )
        disRes$transcriptStart <-
            (disRes$orf_aa_start  * 3 - 2) +
            orfStartDF[
                match(
                    x = disRes$isoform_id,
                    table = orfStartDF$isoform_id
                ),
                2] - 1
        disRes$transcriptEnd <-
            (disRes$orf_aa_end * 3) +
            orfStartDF[
                match(
                    x = disRes$isoform_id,
                    table = orfStartDF$isoform_id
                ),
                2] - 1

        ### convert from transcript to genomic coordinats
        # extract exon data
        myExons <-
            as.data.frame(switchAnalyzeRlist$exons[which(
                switchAnalyzeRlist$exons$isoform_id %in% disRes$isoform_id
            ), ])
        myExonsSplit <- split(myExons, f = myExons$isoform_id)

        # loop over the individual transcripts and extract the genomic coordiants of the domain and also for the active residues (takes 2 min for 17000 rows)
        disResDf <-
            plyr::ddply(
                disRes,
                .progress = progressBar,
                .variables = 'isoform_id',
                .fun = function(aDF) {
                    # aDF <- disRes[which(disRes$isoform_id == 'uc001isa.1'),]

                    transcriptId <- aDF$isoform_id[1]
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

                        orfPosList[[as.character(j)]] <- domainInfo

                    }
                    orfPosDf <- do.call(rbind, orfPosList)

                    return(cbind(aDF, orfPosDf))
                }
            )

        colnames(disResDf) <- gsub(
            'pfam',
            'idr',
            colnames(disResDf)
        )

        disResDf$length <- NULL
    }

    ### Add analysis to switchAnalyzeRlist
    if (TRUE) {
        # sort
        disResDf <-
            disResDf[order(
                disResDf$isoform_id,
                disResDf$transcriptStart
            ), ]

        #disResDf$idrStarExon <- NULL
        #disResDf$idrEndExon <- NULL

        # add the results to the switchAnalyzeRlist object
        switchAnalyzeRlist$idrAnalysis <- disResDf

        # add indication to transcriptDf
        switchAnalyzeRlist$isoformFeatures$idr_identified <- 'no'
        switchAnalyzeRlist$isoformFeatures$idr_identified[which(
            is.na(switchAnalyzeRlist$isoformFeatures$PTC)
        )] <- NA # sets NA for those not analyzed

        switchAnalyzeRlist$isoformFeatures$idr_identified[which(
            switchAnalyzeRlist$isoformFeatures$isoform_id %in%
                disResDf$isoform_id
        )] <- 'yes'
    }

    ### Repport result
    if(TRUE) {
        n <- length(unique(disResDf$isoform_id))
        p <-
            round(n / length(unique(
                switchAnalyzeRlist$isoformFeatures$isoform_id
            )) * 100, digits = 2)

        if (!quiet) {
            message(paste(
                'Added IDR information to ',
                n,
                ' (',
                p,
                '%) transcripts',
                sep = ''
            ))
        }
    }

    return(switchAnalyzeRlist)
}

analyzeDeepLoc2 <- function(
    ### Core arguments
    switchAnalyzeRlist,
    pathToDeepLoc2resultFile,
    
    ### Advanced arguments
    enforceProbabilityCutoff = TRUE,
    probabilityCutoff = NULL,
    ignoreAfterBar = TRUE,
    ignoreAfterSpace = TRUE,
    ignoreAfterPeriod = FALSE,
    quiet = FALSE
) {
    ### Test input
    if(TRUE) {
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(
                'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
            )
        }
        # file
        if (class(pathToDeepLoc2resultFile) != 'character') {
            stop(
                'The \'pathToDeepLoc2resultFile\' argument must be a string pointing to the DeepLoc2 result file(s)'
            )
        }
        if ( ! all(sapply(pathToDeepLoc2resultFile, file.exists)) ) {
            stop('The file(s) \'pathToDeepLoc2resultFile\' points to does not exist')
        }
    }

    ### Import result data
    if(TRUE) {
        ### Test whether it is a tsv
        suppressMessages(
            tmp <- readr::read_tsv( file = pathToDeepLoc2resultFile[1], n_max = 5 )
        )

        ### Read in and test col names
        if(TRUE) {
            ### Tap separated
            if(ncol(tmp) == 13 ) {
                deepLocRes <- do.call(rbind, plyr::llply(
                    pathToDeepLoc2resultFile,
                    .fun = function(
                        aFile
                    ) {
                        suppressMessages(
                            localRes <- readr::read_tsv(aFile)
                        )
                    }
                ))
            }

            ### CSV
            if(ncol(tmp) == 1) {
                deepLocRes <- do.call(rbind, plyr::llply(
                    pathToDeepLoc2resultFile,
                    .fun = function(
                        aFile
                    ) {
                        suppressMessages(
                            localRes <- readr::read_csv(aFile)
                        )
                    }
                ))

            }

            if(! ncol(tmp) %in% c(1, 13) ) {
                ### File format not recogniced
                stop('The file does not appear to be a DeepLoc result file. Double check and try agian...')
            }

            ### Massage result
            if(TRUE) {
                ### remove potential dupliations
                deepLocRes <- dplyr::distinct(deepLocRes)

                ### Test columns
                expectedCols <- c(
                    "Protein_ID",
                    "Localizations",
                    "Signals",
                    "Cytoplasm",
                    "Nucleus",
                    "Extracellular",
                    "Cell membrane",
                    "Mitochondrion",
                    "Plastid",
                    "Endoplasmic reticulum",
                    "Lysosome/Vacuole",
                    "Golgi apparatus",
                    "Peroxisome"
                )

                if( ! all( colnames(deepLocRes) == expectedCols ) ) {
                    stop('The file(s) pointed to by \'pathToDeepLoc2resultFile\' appears not to be DeepLoc Result file')
                }

            }

        }

        ### Test what file contain
        if (nrow(deepLocRes) == 0) {
            stop('The file(s) pointed to by \'pathToDeepLoc2resultFile\' is empty')
        }

        ### Rename
        colnames(deepLocRes)[1] <- 'isoform_id'
        colnames(deepLocRes) <- stringr::str_replace_all( colnames(deepLocRes) , ' ','_')
        colnames(deepLocRes) <- stringr::str_replace_all( colnames(deepLocRes) , '/','_')

        ### Fix names
        if( ignoreAfterBar | ignoreAfterSpace | ignoreAfterPeriod){
          deepLocRes$isoform_id <- fixNames(
            nameVec = deepLocRes$isoform_id,
            ignoreAfterBar = ignoreAfterBar,
            ignoreAfterSpace = ignoreAfterSpace,
            ignoreAfterPeriod = ignoreAfterPeriod
          )
        }
        
        ### Subset to analyzed files
        deepLocRes <- deepLocRes[which(
            deepLocRes$isoform_id %in% switchAnalyzeRlist$orfAnalysis$isoform_id[which(
                !is.na(switchAnalyzeRlist$orfAnalysis$orfTransciptStart)
            )]
        ),]

        if (nrow(deepLocRes) == 0) {
            stop(
                paste(
                    'The file(s) pointed to by \'pathToDeepLoc2resultFile\' does not',
                    'contain results for the isoforms stored in the switchAnalyzeRlist supplied'
                )
            )
        }

        ### Remove singal peptides - use SignalP instead
        deepLocRes$Signals <- NULL

        ### Massage predictions
        deepLocRes$Localizations <- stringr::str_replace_all(deepLocRes$Localizations, ' ','_')
        deepLocRes$Localizations <- stringr::str_replace_all(deepLocRes$Localizations, '/','_')
        deepLocRes$Localizations <- stringr::str_replace_all(deepLocRes$Localizations,'\\|',',')
    }

    ### Filter
    if(TRUE) {
        ### Location specific cutoff
        if(enforceProbabilityCutoff) {
            ### Define location specific probability cutoffs
            probCutoff <- c(
                0.4695,
                0.5275,
                0.6464,
                0.5237,
                0.6373,
                0.6586,
                0.6278,
                0.5648,
                0.6678,
                0.7168
            )
            names(probCutoff) <- c(
                'Cytoplasm',
                'Nucleus',
                'Extracellular',
                'Cell_membrane',
                'Mitochondrion',
                'Plastid',
                'Endoplasmic_reticulum',
                'Lysosome_Vacuole',
                'Golgi_apparatus',
                'Peroxisome'
            )


            colOi <- setdiff(
                colnames(deepLocRes),
                c('isoform_id','Localizations')
            )

            probCutoff <- probCutoff[match(
                colOi, names(probCutoff)
            )]

            ### Identify locations above the specific cutoffs
            deepLocRes$Localizations <- apply(
                deepLocRes[,colOi],
                1,
                function(x) {
                    paste(colOi[which(x >= probCutoff)], collapse = ',')
                }
            )

            ### Set those without location to NA
            deepLocRes$Localizations[which(
                deepLocRes$Localizations == ''
            )] <- NA
        }

        ### User specific cutoff
        if(!is.null(probabilityCutoff)) {
            colOi <- setdiff(
                colnames(deepLocRes),
                c('isoform_id','Localizations')
            )

            ### Identify locations above the specific cutoffs
            deepLocRes$Localizations <- apply(
                deepLocRes[,colOi],
                1,
                function(x) {
                    paste(colOi[which(x >= probabilityCutoff)], collapse = ',')
                }
            )

            ### Set those without location to NA
            deepLocRes$Localizations[which(
                deepLocRes$Localizations == ''
            )] <- NA
        }
    }

    ### Add analysis to switchAnalyzeRlist
    if (TRUE) {
        # add the pfam results to the switchAnalyzeRlist object
        switchAnalyzeRlist$subCellLocationAnalysis <-
            deepLocRes[sort.list(deepLocRes$isoform_id),]

        # add location to isoformFeatures
        switchAnalyzeRlist$isoformFeatures$sub_cell_location <- deepLocRes$Localizations[match(
            switchAnalyzeRlist$isoformFeatures$isoform_id, deepLocRes$isoform_id
        )]
    }


    ### Repport summary
    if(!quiet) {
        nLoc <- sum(!is.na(deepLocRes$Localizations))

        pLoc <-
            round(nLoc / length(unique(
                switchAnalyzeRlist$isoformFeatures$isoform_id
            )) * 100, digits = 2)

        message(paste(
            'Added subcellular information to ', nLoc, ' (', pLoc, '%) transcripts',
            sep = ''
        ))

    }
    return(switchAnalyzeRlist)
}

analyzeDeepTMHMM <- function(
    ### Core arguments
    switchAnalyzeRlist,
    pathToDeepTMHMMresultFile,
    
    ### Advanced arguments
    ignoreAfterBar = TRUE,
    ignoreAfterSpace = TRUE,
    ignoreAfterPeriod = FALSE,
    showProgress = TRUE,
    quiet = FALSE
) {
    ### Test input
    if(TRUE) {
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(
                'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
            )
        }
        if (is.null(switchAnalyzeRlist$orfAnalysis)) {
            stop('ORF needs to be analyzed. Please run \'addORFfromGTF()\' (and if nessesary \'analyzeNovelIsoformORF()\') and try again.')
        }

        # file
        if (class(pathToDeepTMHMMresultFile) != 'character') {
            stop(
                'The \'pathToDeepTMHMMresultFile\' argument must be a string pointing to the DeepTMHMM result file'
            )
        }
        if (! all( file.exists(pathToDeepTMHMMresultFile)) ) {
            stop('(At least on of) the file(s) \'pathToDeepTMHMMresultFile\' points to does not exist')
        }

    }

    if (showProgress & !quiet) {
        progressBar <- 'text'
    } else {
        progressBar <- 'none'
    }

    ### Read result file
    if(TRUE) {
        if (!quiet) {
            message('Step 1 of 2: Reading results into R...')
        }

        ### Read in file
        deepTmRes <- do.call(
            c,
            lapply(pathToDeepTMHMMresultFile, readLines)
        )

        ### Remove unwanted lines
        deepTmRes <- deepTmRes[which(
            ! grepl('^#', deepTmRes)
        )]
        deepTmRes <- deepTmRes[which(
            ! grepl('^/', deepTmRes)
        )]

        ### Read as tsv
        deepTmRes <- read.table(
            text = deepTmRes,
            header = FALSE,
            sep = '\t'
        )

        ### Extract parts of interest
        deepTmRes <- deepTmRes[,1:4]

        colnames(deepTmRes) <- c(
            'isoform_id',
            'region_type',
            'orf_aa_start',
            'orf_aa_end'
        )

        deepTmRes$length <- deepTmRes$orf_aa_end - deepTmRes$orf_aa_start + 1
        
        ### Fix names
        if( ignoreAfterBar | ignoreAfterSpace | ignoreAfterPeriod){
          deepTmRes$isoform_id <- fixNames(
            nameVec = deepTmRes$isoform_id,
            ignoreAfterBar = ignoreAfterBar,
            ignoreAfterSpace = ignoreAfterSpace,
            ignoreAfterPeriod = ignoreAfterPeriod
          )
        }
        
        ### Subset to those in switchList
        deepTmRes <- deepTmRes[which(
            deepTmRes$isoform_id %in%
                switchAnalyzeRlist$orfAnalysis$isoform_id[which(
                    !is.na(switchAnalyzeRlist$orfAnalysis$orfTransciptStart)
                )]
        ), ]
    }
    # deepTmRes

    ### Convert from AA coordinats to transcript and genomic coordinats
    if (TRUE) {
        if (!quiet) {
            message('Step 2 of 2: Converting AA coordinats to transcript and genomic coordinats...')
        }

        ### convert from codons to transcript position
        orfStartDF <-
            unique(
                switchAnalyzeRlist$orfAnalysis[
                    which( !is.na(switchAnalyzeRlist$orfAnalysis$orfTransciptStart)),
                    c('isoform_id', 'orfTransciptStart')
                ]
            )
        deepTmRes$transcriptStart <-
            (deepTmRes$orf_aa_start  * 3 - 2) +
            orfStartDF$orfTransciptStart[
                match(
                    x = deepTmRes$isoform_id,
                    table = orfStartDF$isoform_id
                )] - 1
        deepTmRes$transcriptEnd <-
            (deepTmRes$orf_aa_end * 3) +
            orfStartDF$orfTransciptStart[
                match(
                    x = deepTmRes$isoform_id,
                    table = orfStartDF$isoform_id
                )] - 1

        ### convert from transcript to genomic coordinats
        # extract exon data
        myExons <-
            as.data.frame(switchAnalyzeRlist$exons[which(
                switchAnalyzeRlist$exons$isoform_id %in% deepTmRes$isoform_id
            ), ])
        myExonsSplit <- split(myExons, f = myExons$isoform_id)

        # loop over the individual transcripts and extract the genomic coordiants of the domain and also for the active residues (takes 2 min for 17000 rows)
        deepTmRes <-
            plyr::ddply(
                deepTmRes,
                .progress = progressBar,
                .variables = 'isoform_id',
                .fun = function(aDF) {
                    # aDF <- deepTmRes[which(deepTmRes$isoform_id == 'TCONS_00000045'),]

                    transcriptId <- aDF$isoform_id[1]
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

                        orfPosList[[as.character(j)]] <- domainInfo

                    }
                    orfPosDf <- do.call(rbind, orfPosList)

                    return(cbind(aDF, orfPosDf))
                }
            )

        colnames(deepTmRes) <- gsub(
            'pfam',
            'region',
            colnames(deepTmRes)
        )
    }

    ### Add analysis to switchAnalyzeRlist
    if (TRUE) {
        # sort
        deepTmRes <-
            deepTmRes[order(
                deepTmRes$isoform_id,
                deepTmRes$transcriptStart
            ), ]

        # add the results to the switchAnalyzeRlist object
        switchAnalyzeRlist$topologyAnalysis <- deepTmRes

    }

    ### Repport result
    if(TRUE) {
        n <- length(unique(deepTmRes$isoform_id))
        p <-
            round(n / length(unique(
                switchAnalyzeRlist$isoformFeatures$isoform_id
            )) * 100, digits = 2)

        if (!quiet) {
            message(paste(
                'Added topology information to ',
                n,
                ' transcripts (',
                p,
                '%).',
                sep = ''
            ))
        }
    }

    return(switchAnalyzeRlist)
}
