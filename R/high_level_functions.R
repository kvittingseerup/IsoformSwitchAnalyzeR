isoformSwitchAnalysisPart1 <- function(
    input,
    alpha=0.05,
    dIFcutoff = 0.1,
    calibratePvalues=TRUE,
    orfMethod='longest',
    genomeObject,
    cds=NULL,
    pathToOutput=getwd(),
    outputSequences=TRUE,
    overwritePvalues=TRUE,
    overwriteORF=FALSE,
    quiet=FALSE
) {
    nrAnalysis <- 4

    ### Identfy input type
    if(TRUE) {
        ### Import cuffdiff data
        inputClass <- class(input)
        if(inputClass == 'switchAnalyzeRlist') {
            if( ! quiet) { message(paste('Step 1 of', nrAnalysis, ': Massaging data...', sep=' ')) }

            switchAnalyzeRlist <- input
        } else if (inputClass == 'CuffSet'){
            if( ! quiet) { message(paste('Step 1 of', nrAnalysis, ': Importing and massaging Cufflinks/Cuffdiff data...', sep=' ')) }

            switchAnalyzeRlist <- importCufflinksData(cuffDB = input, quiet=TRUE)
        } else if (inputClass == 'ballgown'){
            if( ! quiet) { message(paste('Step 1 of', nrAnalysis, ': Importing and massaging data from ballgown...', sep=' ')) }

            switchAnalyzeRlist <- importBallgownData(ballgownObject = input, quiet = TRUE)
        } else {
            stop('Input is not recogniced - see ?isoformSwitchAnalysisPart1 for more details')
        }
    }

    ### Run preFilter
    switchAnalyzeRlist <- preFilter(switchAnalyzeRlist=switchAnalyzeRlist, removeSingleIsoformGenes = TRUE, quiet=TRUE)

    ### Test isoform switches
    noTestPerformed <- all(is.na( switchAnalyzeRlist$isoformFeatures$gene_switch_q_value ))
    if( noTestPerformed | overwritePvalues ) {
        if( ! quiet) {message(paste('Step 2 of', nrAnalysis, ': Detecting isoform switches...', sep=' ')) }
        switchAnalyzeRlist <- isoformSwitchTest(switchAnalyzeRlist = switchAnalyzeRlist, reduceToSwitchingGenes = TRUE, alpha = alpha, dIFcutoff=dIFcutoff, calibratePvalues = calibratePvalues, quiet=TRUE)
    } else {
        if( ! quiet) {message(paste('Step 2 of', nrAnalysis, ': No isoform switch was performed...', sep=' ')) }
    }

    ### Predict ORF
    if( ! quiet) { message(paste('Step 3 of', nrAnalysis, ': Predicting open reading frames', sep=' ')) }
    if( all(names(switchAnalyzeRlist) !='orfAnalysis') | overwriteORF){
        switchAnalyzeRlist <- analyzeORF(switchAnalyzeRlist = switchAnalyzeRlist, genomeObject = genomeObject, cds = cds, orfMethod = orfMethod, quiet=TRUE)
    }

    ### Extract and write sequences
    if( ! quiet) { message(paste('Step 4 of', nrAnalysis, ': Extracting (and outputting) sequences', sep=' ')) }
    switchAnalyzeRlist <- extractSequence(
        switchAnalyzeRlist = switchAnalyzeRlist,
        genomeObject = genomeObject,
        onlySwitchingGenes = TRUE,
        extractNTseq = TRUE,
        extractAAseq = TRUE,
        addToSwitchAnalyzeRlist = TRUE,
        writeToFile = outputSequences,
        pathToOutput = pathToOutput,
        filterAALength=TRUE,
        quiet=TRUE
    )

    ### Print summary
    if( ! quiet) {
        message('\nThe number of isoform swithces found were:')
        print(extractSwitchSummary(switchAnalyzeRlist, alpha = alpha, dIFcutoff=dIFcutoff))
        message(
            paste(
                'The nucleotide and amino acid sequences of these isoforms have been outputted',
                'to the supplied directory enabling external analysis of protein domians (Pfam), coding potential (CPAT) or signal peptides (SignalIP).',
                '\nSee ?analyzeCPAT, ?analyzePFAM or ?analyzeSignalIP (under details) for suggested ways of running these three tools.',
                sep=' '
            )
        )
    }
    return(switchAnalyzeRlist)
}

isoformSwitchAnalysisPart2 <- function(
    switchAnalyzeRlist,
    alpha=0.05,
    dIFcutoff = 0.1,
    n=NA,
    codingCutoff,
    removeNoncodinORFs,
    pathToCPATresultFile=NULL,
    pathToPFAMresultFile=NULL,
    pathToSignalPresultFile=NULL,
    consequencesToAnalyze=c('intron_retention','coding_potential','ORF_seq_similarity','NMD_status','domains_identified','signal_peptide_identified'),
    pathToOutput=getwd(),
    fileType='pdf',
    asFractionTotal=FALSE,
    outputPlots=TRUE,
    quiet=FALSE
) {
    ### Test input
    if(TRUE) {
        if(!is.null(pathToCPATresultFile)) {
            if(is.na(codingCutoff)) {
                stop(
                    paste(
                        'A cutoff must be supplied to codingCutoff if a CPAT analysis are added to pathToCPATresultFile.',
                        'The cutoff is dependent on spieces analyzed. The suggested cutoffs from the CPAT article (see refrences) is 0.364 for human and 0.44 for mouse.',
                        'see ?analyzeCPAT for more information.',
                        sep=' '
                        )
                )
            }
        }
    }

    nrAnalysis <- 2 + as.integer( any(c('all','intron_retention') %in% consequencesToAnalyze) ) + 2*as.integer(outputPlots)
    analysisDone <- 1

    ### Add annoation
    if( ! quiet) { message(paste('Step', analysisDone,  'of', nrAnalysis, ': Importing external sequence analysis...', sep=' ')) }

    if(!is.null(pathToCPATresultFile)) {
        switchAnalyzeRlist <- analyzeCPAT(switchAnalyzeRlist = switchAnalyzeRlist, pathToCPATresultFile = pathToCPATresultFile, codingCutoff = codingCutoff, removeNoncodinORFs=removeNoncodinORFs, quiet=TRUE)
    }
    if(!is.null(pathToPFAMresultFile)) {
        switchAnalyzeRlist <- analyzePFAM(switchAnalyzeRlist = switchAnalyzeRlist, pathToPFAMresultFile = pathToPFAMresultFile, quiet=TRUE)
    }
    if(!is.null(pathToSignalPresultFile)) {
        switchAnalyzeRlist <- analyzeSignalP(switchAnalyzeRlist = switchAnalyzeRlist, pathToSignalPresultFile = pathToSignalPresultFile, quiet=TRUE)
    }
    analysisDone <- analysisDone +1

    ### Predict intron retentions
    if( any(c('all','intron_retention') %in% consequencesToAnalyze)) {
        if( ! quiet) { message(paste('Step', analysisDone,  'of', nrAnalysis, ': Analyzing intron retentions...', sep=' ')) }

        switchAnalyzeRlist <- analyzeIntronRetention(switchAnalyzeRlist = switchAnalyzeRlist, onlySwitchingGenes = TRUE, alpha = alpha, dIFcutoff=dIFcutoff, quiet=TRUE)
        analysisDone <- analysisDone +1
    }

    ### Predict functional consequences
    if( ! quiet) { message(paste('Step', analysisDone,  'of', nrAnalysis, ': Prediciting functional consequences...', sep=' ')) }

    switchAnalyzeRlist <- analyzeSwitchConsequences(switchAnalyzeRlist = switchAnalyzeRlist, consequencesToAnalyze = consequencesToAnalyze, alpha = alpha, dIFcutoff = dIFcutoff, quiet=TRUE)
    analysisDone <- analysisDone +1

    ### Make isoform switch plots
    if( outputPlots ) {
        if( ! quiet) { message(paste('Step', analysisDone,  'of', nrAnalysis, ': Making indidual isoform switch plots...', sep=' ')) }

        switchPlotTopSwitches(switchAnalyzeRlist = switchAnalyzeRlist, alpha = alpha, dIFcutoff=dIFcutoff, n = n, pathToOutput = pathToOutput, filterForConsequences=TRUE, splitFunctionalConsequences = FALSE, fileType = fileType, quiet=TRUE)

        analysisDone <- analysisDone +1
    }

    ### Make overall consequences
    if( outputPlots ) {
        if( ! quiet) { message(paste('Step', analysisDone,  'of', nrAnalysis, ': Analyzing combined consequences plot...', sep=' ')) }
        if(fileType == 'pdf') {
            pdf(file = paste(pathToOutput,'/common_switch_consequences.pdf', sep=''), width = 10, height = 7)
            extractConsequenceSummary(switchAnalyzeRlist = switchAnalyzeRlist, asFractionTotal = asFractionTotal, alpha = alpha, dIFcutoff=dIFcutoff, plot = TRUE, returnResult = FALSE)
            dev.off()
        } else {
            png(file = paste(pathToOutput,'/common_switch_consequences.png', sep=''), width = 1000, height = 700)
            extractConsequenceSummary(switchAnalyzeRlist = switchAnalyzeRlist, asFractionTotal = asFractionTotal, alpha = alpha, dIFcutoff=dIFcutoff, plot = TRUE, returnResult = FALSE)
            dev.off()
        }
    }

    ### Print summary
    if( ! quiet) {
        message('\nThe number of isoform swithces with functional consequences identified were:')
        print(extractSwitchSummary(switchAnalyzeRlist, alpha = alpha, dIFcutoff=dIFcutoff, filterForConsequences = TRUE))
        if( outputPlots ) {
            message(
                paste(
                    'The switch analysis plot for each of these, as well as a plot summarizing the functional consequences',
                    'have been outputted to the folder specified by \'pathToOutput\'.',
                    sep=' '
                )
            )
        }
    }

    return(switchAnalyzeRlist)
}

isoformSwitchAnalysisCombined <- function(
    input,
    alpha=0.05,
    dIFcutoff = 0.1,
    calibratePvalues=TRUE,
    n=NA,
    pathToOutput=getwd(),
    overwritePvalues=TRUE,
    overwriteORF=FALSE,
    outputSequences=FALSE,
    genomeObject,
    orfMethod='longest',
    cds=NULL,
    consequencesToAnalyze=c('intron_retention','ORF_seq_similarity','NMD_status'),
    fileType='pdf',
    asFractionTotal=FALSE,
    outputPlots=TRUE,
    quiet=FALSE
) {
    ### Run part 1
    if( ! quiet) { message('\nPART 1: EXTRACTING ISOFORM SWTICH SEQUENCES') }
    switchAnalyzeRlist <- isoformSwitchAnalysisPart1(input = input, alpha = alpha, dIFcutoff=dIFcutoff, calibratePvalues = calibratePvalues, pathToOutput = pathToOutput, genomeObject = genomeObject, orfMethod = orfMethod, cds = cds, outputSequences=outputSequences, overwritePvalues=overwritePvalues, overwriteORF=overwriteORF, quiet=quiet)

    ### Run part 2 without annoation
    if( ! quiet) { message('\nPART 2: PLOTTING ISOFORM SEQUENCES') }
    switchAnalyzeRlist <- isoformSwitchAnalysisPart2(switchAnalyzeRlist = switchAnalyzeRlist, alpha = alpha, dIFcutoff=dIFcutoff, n=n, consequencesToAnalyze=consequencesToAnalyze, pathToOutput=pathToOutput, fileType=fileType, asFractionTotal=asFractionTotal, outputPlots=outputPlots, quiet=quiet)

    return(switchAnalyzeRlist)
}
