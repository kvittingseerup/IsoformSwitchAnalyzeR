extractExpressionMatrix <- function(
    switchAnalyzeRlist,
    feature='isoformUsage',
    addInfo=FALSE,
    na.rm=TRUE
) {
    ### Test input
    if( class(switchAnalyzeRlist) != 'switchAnalyzeRlist' )        { stop('The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\'') }
    if( ! feature %in% c('geneExp','isoformExp','isoformUsage')) { stop('The \'feature\' argument must be either \'geneExp\' or \'isoformExp\' or \'IF\'.')}

    if(feature == 'geneExp') {
        localData1 <- unique( switchAnalyzeRlist$isoformFeatures[,c('gene_id','condition_1','gene_value_1')] )
        localData2 <- unique( switchAnalyzeRlist$isoformFeatures[,c('gene_id','condition_2','gene_value_2')] )

        recast1 <- reshape2::dcast(localData1, gene_id ~ condition_1, value.var='gene_value_1')
        recast2 <- reshape2::dcast(localData2, gene_id ~ condition_2, value.var='gene_value_2')

        notIn1 <- recast2[,c('gene_id',setdiff(colnames(recast2), colnames(recast1)))]

        combinedData <- unique( merge(recast1, notIn1, by='gene_id') )

        if(na.rm) {
            combinedData <- combinedData[ which( apply(combinedData, 1, function(x) { ! any(is.na(x)) } )),]
        }

        # potentially add info
        if(addInfo) {
            # extract info
            geneInfo <- unique(switchAnalyzeRlist$isoformFeatures[, which(
                colnames(switchAnalyzeRlist$isoformFeatures) %in%
                    c('gene_id','gene_name')
            )])
            # merge
            combinedData <- merge(combinedData, geneInfo, by='gene_id')
        } else {
            rownames(combinedData) <- combinedData$gene_id
            combinedData$gene_id <- NULL
        }

    } else if( feature == 'isoformExp') {
        localData1 <- unique( switchAnalyzeRlist$isoformFeatures[,c('isoform_id','condition_1','iso_value_1')] )
        localData2 <- unique( switchAnalyzeRlist$isoformFeatures[,c('isoform_id','condition_2','iso_value_2')] )

        recast1 <- reshape2::dcast(localData1, isoform_id ~ condition_1, value.var='iso_value_1')
        recast2 <- reshape2::dcast(localData2, isoform_id ~ condition_2, value.var='iso_value_2')

        notIn1 <- recast2[,c('isoform_id',setdiff(colnames(recast2), colnames(recast1)))]

        combinedData <- unique( merge(recast1, notIn1, by='isoform_id') )

        if(na.rm) {
            combinedData <- combinedData[ which( apply(combinedData, 1, function(x) { ! any(is.na(x)) } )),]
        }

        # potentially add info
        if(addInfo) {
            # extract iso info
            isoInfo <- unique(switchAnalyzeRlist$isoformFeatures[, which(
                colnames(switchAnalyzeRlist$isoformFeatures) %in%
                    c('isoform_id','gene_id','gene_name','nearest_ref_id','class_code','length','IR','signal_peptide_identified',
                      'codingPotentialValue','codingPotential','domain_identified','subcellularOrign')
            )])
            # extract ORF
            orfInfo <- unique(switchAnalyzeRlist$orfAnalysis[, which(
                colnames(switchAnalyzeRlist$orfAnalysis) %in%
                    c('isoform_id','orfTransciptStart','orfTransciptEnd','orfTransciptLength','orfStartGenomic','orfEndGenomic','PTC')
            )])
            # merge
            isoInfo <- merge(isoInfo, orfInfo, by='isoform_id')

            # merge
            combinedData <- merge(combinedData, isoInfo, by='isoform_id')
        } else {
            rownames(combinedData) <- combinedData$isoform_id
            combinedData$isoform_id <- NULL
        }
    } else {
        localData1 <- unique( switchAnalyzeRlist$isoformFeatures[,c('isoform_id','condition_1','IF1')] )
        localData2 <- unique( switchAnalyzeRlist$isoformFeatures[,c('isoform_id','condition_2','IF2')] )

        recast1 <- reshape2::dcast(localData1, isoform_id ~ condition_1, value.var='IF1')
        recast2 <- reshape2::dcast(localData2, isoform_id ~ condition_2, value.var='IF2')

        notIn1 <- recast2[,c('isoform_id',setdiff(colnames(recast2), colnames(recast1)))]

        combinedData <- unique( merge(recast1, notIn1, by='isoform_id') )

        if(na.rm) {
            combinedData <- combinedData[ which( apply(combinedData, 1, function(x) { ! any(is.na(x)) } )),]
        }

        # potentially add info
        if(addInfo) {
            # extract iso info
            isoInfo <- unique(switchAnalyzeRlist$isoformFeatures[, which(
                colnames(switchAnalyzeRlist$isoformFeatures) %in%
                    c('isoform_id','gene_id','gene_name','nearest_ref_id','class_code','length','IR','signal_peptide_identified',
                      'codingPotentialValue','codingPotential','domain_identified','subcellularOrign')
            )])
            # extract ORF
            orfInfo <- unique(switchAnalyzeRlist$orfAnalysis[, which(
                colnames(switchAnalyzeRlist$orfAnalysis) %in%
                    c('isoform_id','orfTransciptStart','orfTransciptEnd','orfTransciptLength','orfStartGenomic','orfEndGenomic','PTC')
            )])
            # merge
            isoInfo <- merge(isoInfo, orfInfo, by='isoform_id')

            # merge
            combinedData <- merge(combinedData, isoInfo, by='isoform_id')
        } else {
            rownames(combinedData) <- combinedData$isoform_id
            combinedData$isoform_id <- NULL
        }
    }

    return(combinedData)
}

prepareCuffExample <- function(){
    dir <- tempdir()
    extdata <- system.file("extdata", package="cummeRbund")
    file.copy(file.path(extdata, dir(extdata)), dir)

    cuffDB <- cummeRbund::readCufflinks(dir=dir, gtfFile=system.file("extdata/chr1_snippet.gtf", package="cummeRbund"), genome="hg19", rebuild=TRUE)
    return(cuffDB)
}
