makeMinimumSwitchList <- function(
    orgSwitchList,
    isoformsToKeep
){
    ### Subset to wanted isoforms
    orgSwitchList <- subset(orgSwitchList, orgSwitchList$isoformFeatures$isoform_id %in% isoformsToKeep)

    ### remove non-needed entries
    orgSwitchList$isoformSwitchAnalysis <- NULL
    orgSwitchList$switchConsequence <- NULL

    ### Reduce size of isoformFeatures
    colsToAnnulate <- c(
        'gene_name','gene_value_1','gene_value_2',
        'gene_stderr_1','gene_stderr_2','gene_log2_fold_change','gene_q_value',
        'iso_value_1','iso_value_2','iso_stderr_1','iso_stderr_2',
        'iso_log2_fold_change','iso_q_value','IF1','IF2','dIF',
        'isoform_switch_q_value')
    orgSwitchList$isoformFeatures <- orgSwitchList$isoformFeatures[, which( ! colnames(orgSwitchList$isoformFeatures) %in% colsToAnnulate)]
    orgSwitchList$isoformFeatures$condition_1 <- 1
    orgSwitchList$isoformFeatures$condition_2 <- 2
    orgSwitchList$isoformFeatures$gene_switch_q_value <- 1
    orgSwitchList$isoformFeatures$isoform_switch_q_value <- 1

    orgSwitchList$isoformFeatures <- unique(orgSwitchList$isoformFeatures)

    return(orgSwitchList)
}

analyzeIntronRetention <- function(
    switchAnalyzeRlist,
    onlySwitchingGenes=TRUE,
    alpha=0.05,
    dIFcutoff = 0.1,
    showProgress=TRUE,
    quiet=FALSE
) {
    ### Test input
    if(TRUE) {
        if( class(switchAnalyzeRlist) != 'switchAnalyzeRlist' ) { stop('The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\'') }
        if( !is.logical(onlySwitchingGenes)) {stop('The onlySwitchingGenes argument must be either TRUE or FALSE')}
        if( onlySwitchingGenes ) {
            if( !any( !is.na(switchAnalyzeRlist$isoformFeatures$gene_switch_q_value)) ) {
                stop('The analsis of isoform switching must be performed before functional consequences can be analyzed. Please run ?detectIsoformSwitching and try again.')
            }
        }
        if( alpha < 0 | alpha > 1 ) { stop('The alpha parameter must be between 0 and 1 ([0,1]).')}
        if( alpha > 0.05) {
            warning('Most journals and scientists consider an alpha larger than 0.05 untrustworthy. We therefore recommend using alpha values smaller than or queal to 0.05')
        }
        if( dIFcutoff < 0 | dIFcutoff > 1) {
            stop('The dIFcutoff must be in the interval [0,1].')
        }
    }

    ### Exract data
    if(TRUE) {
        message('Step 1 of 4: Massaging data...')
        localData <- unique(
            makeMinimumSwitchList(switchAnalyzeRlist, switchAnalyzeRlist$isoformFeatures$isoform_id)$isoformFeatures[,
                c('isoform_id','gene_id','condition_1','condition_2')
            ]
        )
        localData$iso_value_1 <- NA
        localData$iso_value_2 <- NA
        localData$iso_q_value <- NA
        localData$gene_value_1 <- NA
        localData$gene_value_2 <- NA

        colnames(localData)[match( c('condition_1','condition_2') , colnames(localData))] <- c('sample_1','sample_2')

        if( onlySwitchingGenes ) {
            isoResTest <- any( ! is.na(switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value) )
            if( isoResTest ) {
                switchingGenes <- unique( switchAnalyzeRlist$isoformFeatures$gene_id [ which(
                    switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value  < alpha     &
                        abs( switchAnalyzeRlist$isoformFeatures$dIF ) > dIFcutoff
                ) ] )
            } else {
                switchingGenes <- unique( switchAnalyzeRlist$isoformFeatures$gene_id [ which(
                    switchAnalyzeRlist$isoformFeatures$gene_switch_q_value  < alpha     &
                        abs( switchAnalyzeRlist$isoformFeatures$dIF ) > dIFcutoff
                ) ] )
            }
            if(length(switchingGenes) == 0) {stop('No switching genes were found. Pleasae turn off \'onlySwitchingGenes\' and try again.')}

            localData <- localData[which(localData$gene_id %in% switchingGenes),]

        }

        correspondingExons <- switchAnalyzeRlist$exons[which( switchAnalyzeRlist$exons$isoform_id %in% localData$isoform_id ), ]
        colnames(mcols(correspondingExons)) <- paste('spliceR.',colnames(mcols(correspondingExons)), sep='')
    }

    ### Make temporary spliceRList
    if(TRUE) {
        message('Step 2 of 4: Creating spliceRList...')
        # make GRange
        n <- nrow(localData)
        transcriptFeatureGR <- GRanges(
            seqnames=rep('NA',n),
            ranges=IRanges(
                start=rep(0, n),
                end=rep(1,n)
            ),
            spliceR=localData		#add spliceR to metadata colnames
        )

        # Make local spliceRList
        localSpliceRList <- new("SpliceRList",
                                list(
                                    transcript_features=transcriptFeatureGR,
                                    exon_features=correspondingExons,
                                    assembly_id='NA',
                                    source_id='cufflinks',
                                    conditions=switchAnalyzeRlist$conditions$condition,
                                    transcripts_plot=NULL,
                                    filter_params=NULL
                                )
        )

    }

    ### Run spliceR::spliceR
    message('Step 3 of 4: Running spliceR...')
    if(quiet) {
        suppressMessages(
            annotationResult <- spliceR::spliceR(
                transcriptData = localSpliceRList,
                compareTo = 'preTranscript',
                filters = 'none',
                useProgressBar = showProgress & !quiet
            )
        )
    } else {
        annotationResult <- spliceR::spliceR(
            transcriptData = localSpliceRList,
            compareTo = 'preTranscript',
            filters = 'none',
            useProgressBar = showProgress & !quiet
        )
    }

    message('Step 4 of 4: Preparing output...')
    localIRresults <- unique( as.data.frame(mcols(annotationResult$transcript_features)[,c('spliceR.isoform_id','spliceR.ISI','spliceR.ISI.start','spliceR.ISI.end')]) )
    colnames(localIRresults) <- gsub('\\.','_', gsub('spliceR.', '', colnames(localIRresults)))
    colnames(localIRresults) <- gsub('ISI', 'IR', colnames(localIRresults))
    colnames(localIRresults)[3:4] <- gsub('IR', 'IR_genomic', colnames(localIRresults)[3:4])

    ### Transfer result to switchAnalyzeRlist list
    switchAnalyzeRlist$isoformFeatures$IR <- localIRresults$IR[
        match( switchAnalyzeRlist$isoformFeatures$isoform_id , localIRresults$isoform_id) ]

    switchAnalyzeRlist$intronRetentionAnalysis <- localIRresults


    return(switchAnalyzeRlist)
}
