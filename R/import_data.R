### Helper functions
allPairwiseFeatures <- function(aNameVec1, forceNonOverlap = FALSE) {
    aNameVec1 <- unique(as.character(aNameVec1))

    # vectors to store result
    var1 <- character()
    var2 <- character()

    ### Get comparisons
    count <- 2
    n <- length(aNameVec1)
    for (i in 1:(n - 1)) {
        for (j in count:n) {
            var1 <- c(var1, aNameVec1[i])
            var2 <- c(var2, aNameVec1[j])
        }
        count <- count + 1
    }

    myCombinations <-
        data.frame(var1 = var1,
                   var2 = var2,
                   stringsAsFactors = FALSE)

    return(myCombinations)
}
jaccardDistance <- function(x, y) {
    length(intersect(x, y)) / length(union(x, y))
}

### Acutal import functions
importCufflinksCummeRbund <- function(
    cuffDB,
    fixCufflinksAnnotationProblem = TRUE,
    addCufflinksSwichTest = TRUE,
    quiet = FALSE
) {
    ### Test input
    if (TRUE) {
        if (class(cuffDB)      != 'CuffSet') {
            stop(paste(
                'The object supplied to \'cuffDB\' must be a \'CuffSet\'.',
                'See ?cummeRbund::readCufflinks for more information'
            ))
        }
        if (!is.logical(fixCufflinksAnnotationProblem)) {
            stop(
            "The \'fixCufflinksAnnotationProblem\' argument must be a logic."
            )
        }

    }
    if (!quiet) {
        message("Reading cuffDB, isoforms...")
    }

    # Make gene and isoform pointers
    cuffGenes 			<- cummeRbund::genes(cuffDB)
    cuffIsoforms 		<- cummeRbund::isoforms(cuffDB)

    ### Extract gene diff analysis
    geneDiffanalysis     		<-
        data.frame(cummeRbund::diffData(cuffGenes), stringsAsFactors = FALSE)
    geneDiffanalysis <-
        geneDiffanalysis[which(colnames(geneDiffanalysis) != 'test_stat')]
    colnames(geneDiffanalysis) 	<-
        c('gene_id', 'sample_1', 'sample_2' ,
          unlist(lapply(
              colnames(geneDiffanalysis[, -1:-3]),
              function(x) {
                  paste("gene_", x, sep = "")
              }))
        ) # add gene to the colnames so they can be destinquished from the isoform diff data

    ### Extract isoform (and thereby also gene) annotation
    isoformAnnotation <- data.frame(
        cummeRbund::annotation(cuffIsoforms), stringsAsFactors = FALSE)

    # Unique + removeal of colums is nessesary because of the new cummeRbund devel (where these colums are included and exon number thereby makes duplication of rows)
    isoformAnnotation <-
        unique(isoformAnnotation[, na.omit(match(
            c(
                'isoform_id',
                'gene_id',
                'gene_short_name',
                'nearest_ref_id',
                'class_code',
                'TSS_group_id',
                'CDS_id',
                'length',
                'locus'
            ),
            colnames(isoformAnnotation)
        ))])

    ### Extract isoform diff analysis
    isoformDiffanalysis <- data.frame(
        cummeRbund::diffData(cuffIsoforms),stringsAsFactors = FALSE)

    isoformDiffanalysis <-
        isoformDiffanalysis[, which(
            colnames(isoformDiffanalysis) != 'test_stat'
        )]
    colnames(isoformDiffanalysis) 	<-
        c('isoform_id', 'sample_1', 'sample_2' ,
          unlist(lapply(colnames(isoformDiffanalysis[, -1:-3]), function(x)
         {paste("iso_", x, sep = "")} ))) # add gene to the colnames so they can be destinquished from the gene diff data

    ### Extract and add isoform stderr
    # exract
    isoStderr <- cummeRbund::fpkm(cuffIsoforms)[, c(
            'isoform_id', 'sample_name', 'stdev'
        )]
    colnames(isoStderr)[which(colnames(isoStderr) == 'stdev')] <- 'iso_stderr'

    # add
    isoformDiffanalysis <-
        merge(
            isoformDiffanalysis,
            isoStderr,
            by.x = c('isoform_id', 'sample_2'),
            by.y = c('isoform_id', 'sample_name')
        )
    colnames(isoformDiffanalysis)[which(
        grepl('iso_stderr', colnames(isoformDiffanalysis)))] <- 'iso_stderr_2'

    isoformDiffanalysis <-
        merge(
            isoformDiffanalysis,
            isoStderr,
            by.x = c('isoform_id', 'sample_1'),
            by.y = c('isoform_id', 'sample_name')
        )
    colnames(isoformDiffanalysis)[ which(
        grepl('iso_stderr$', colnames(isoformDiffanalysis), perl = TRUE)
    )] <- c('iso_stderr_1')

    # reorder
    isoformDiffanalysis <-
        isoformDiffanalysis[, c(
            'isoform_id',
            'sample_1',
            'sample_2',
            'iso_status',
            'iso_value_1',
            'iso_value_2',
            'iso_stderr_1',
            'iso_stderr_2',
            'iso_log2_fold_change',
            'iso_p_value',
            'iso_q_value',
            'iso_significant'
        )]

    ### Extract and add gene stderr
    # extract
    geneStderr <-
        fpkm(cuffGenes)[, c(c('gene_id', 'sample_name', 'stdev'))]
    colnames(geneStderr)[which(colnames(geneStderr) == 'stdev')] <-
        'gene_stderr'

    # Add
    geneDiffanalysis <-
        merge(
            geneDiffanalysis,
            geneStderr,
            by.x = c('gene_id', 'sample_2'),
            by.y = c('gene_id', 'sample_name')
        )
    colnames(geneDiffanalysis)[which(
        grepl('gene_stderr', colnames(geneDiffanalysis)))] <- 'gene_stderr_2'

    geneDiffanalysis <-
        merge(
            geneDiffanalysis,
            geneStderr,
            by.x = c('gene_id', 'sample_1'),
            by.y = c('gene_id', 'sample_name')
        )
    colnames(geneDiffanalysis)[which(
        grepl('gene_stderr$', colnames(geneDiffanalysis), perl = TRUE)
    )] <- c('gene_stderr_1')

    # reorder
    geneDiffanalysis <-
        geneDiffanalysis[, c(
            'gene_id',
            'sample_1',
            'sample_2',
            'gene_status',
            'gene_value_1',
            'gene_value_2',
            'gene_stderr_1',
            'gene_stderr_2',
            'gene_log2_fold_change',
            'gene_p_value',
            'gene_q_value',
            'gene_significant'
        )]

    ## Merge isoform annotation and gene diff analysis
    isoformData <-
        merge(isoformAnnotation, geneDiffanalysis, by = 'gene_id')
    ## Merge with isoform diff analysis
    isoformData <-
        merge(isoformData,
              isoformDiffanalysis,
              by = c('isoform_id', 'sample_1', 'sample_2'))

    colnames(isoformData)[which(
        colnames(isoformData) == 'gene_short_name'
    )] <- 'gene_name'


    ### Extract exon info
    if (!quiet) {
        message("Reading cuffDB, exons...")
    }
    isoformFeatureQuery		<-
        paste(
            "SELECT y.* FROM features y JOIN genes x on y.gene_id = x.gene_id ",
              sep = ""
        )
    isoformFeatures			<-
        data.frame(dbGetQuery(cuffDB@DB, isoformFeatureQuery),
                   stringsAsFactors = FALSE)
    isoformFeatures         <-
        isoformFeatures[which(isoformFeatures$strand %in% c('+', '-')), ]
    isoformFeatures         <-
        isoformFeatures [which(isoformFeatures$type == 'exon'), ]

    if (nrow(isoformFeatures) == 0) {
        stop(
            "No exon information extracted - this is moste likely because the cuffDB was not build with a GTF file."
        )
    }
    # Another faster way of doing it would be to use the colums removed from isoformAnnotation, but this approach is not backwards compatible
    # isoformFeatures <- isoformAnnotation[,c("seqnames", "start", "end", "strand", "isoform_id", "gene_id")]

    ### Extract isoform replicate expression
    isoRepExp <- repCountMatrix(cuffIsoforms)
    isoRepExp$isoform_id <- rownames(isoRepExp)
    rownames(isoRepExp) <- NULL
    isoRepExp <-
        isoRepExp[, c(ncol(isoRepExp), 1:(ncol(isoRepExp) - 1))]

    ### Make sure both all data.frames only contain isoforms also found in the other data.frame
    if (TRUE) {
        myUnion     <-
            unique(
                c(
                    isoformData$isoform_id,
                    isoformFeatures$isoform_id,
                    isoRepExp$isoform_id
                )
            )
        myIntersect <- intersect(
            intersect(isoformData$isoform_id, isoformFeatures$isoform_id),
            isoRepExp$isoform_id
        )

        # If there are descripencies
        if (length(myUnion) != length(myIntersect)) {
            isoformData      <- isoformData[which(
                isoformData$isoform_id %in% myIntersect
            ),]
            isoformFeatures  <- isoformFeatures[which(
                isoformFeatures$isoform_id %in% myIntersect
            ), ]
            isoRepExp <- isoRepExp[which(
                isoRepExp$isoform_id %in% myIntersect
            ), ]

            if (!quiet) {
                message(
                    paste(
                        'There were discrepencies between the GTF and',
                        'the expression analysis files. To solve this',
                        abs(length(myUnion) - length(myIntersect)) ,
                        'transcripts were removed.',
                        sep = ' '
                    )
                )
            }
        }

    }

    ### Fix to correct for Cufflinks annotation problem where cufflinks
    # assignes transcripts from several annotated genes to 1 cuffgene
    if (fixCufflinksAnnotationProblem) {
        if (!quiet) {
            message("Analyzing cufflinks annotation problem...")
        }

        geneName <- unique(isoformData[, c('gene_id', 'gene_name')])
        geneNameSplit <-
            split(geneName$gene_name , f = geneName$gene_id)
        # remove all unique
        geneNameSplit <-
            geneNameSplit[which(sapply(geneNameSplit, function(x)
                length(unique(x))) > 1)]

        if (length(geneNameSplit) > 0) {
            # if there are any problems
            if (!quiet) {
                message("Fixing cufflinks annotation problem...")
            }
            #get indexes of those affected
            geneNameIndexesData     <-
                which(isoformData$gene_id %in% names(geneNameSplit))
            geneNameIndexesFeatures <-
                which(isoformFeatures$gene_id %in% names(geneNameSplit))
            # combine names of cuffgenes and
            isoformData$gene_id[geneNameIndexesData]            <-
                paste(isoformData$gene_id[geneNameIndexesData]        ,
                      isoformData$gene_name[geneNameIndexesData],
                      sep = ':')
            isoformFeatures$gene_id[geneNameIndexesFeatures]    <-
                paste(isoformFeatures$gene_id[geneNameIndexesFeatures],
                      isoformFeatures$gene_name[geneNameIndexesFeatures],
                      sep = ':')

            ## Correct gene expression levels and differntial analysis
            problematicGenes <-
                isoformData[geneNameIndexesData, c(
                    'isoform_id',
                    'gene_id',
                    'sample_1',
                    'sample_2',
                    'gene_value_1',
                    'gene_value_2',
                    'gene_stderr_1',
                    'gene_stderr_2',
                    'gene_log2_fold_change',
                    'gene_p_value',
                    'gene_q_value',
                    'gene_significant',
                    'iso_status',
                    'iso_value_1',
                    'iso_value_2'
                )]
            problematicGenesSplit <-
                split(problematicGenes, f = problematicGenes[, c('gene_id', 'sample_1', 'sample_2')], drop =
                          TRUE)

            correctedGenes <-
                ldply(
                    problematicGenesSplit,
                    .fun = function(df) {
                        # df <- problematicGenesSplit[[1]]
                        df$gene_value_1 <- sum(df$iso_value_1)
                        df$gene_value_2 <- sum(df$iso_value_2)
                        df$gene_stderr_1 <- NA
                        df$gene_stderr_2 <- NA
                        df$gene_log2_fold_change <-
                            log2((df$gene_value_2[2] + 0.0001) / (df$gene_value_1[1] + 0.0001))
                        df$gene_p_value <- 1
                        df$gene_q_value <- 1
                        df$gene_significant <- 'no'
                        df$iso_status <- 'NOTEST'
                        return(df)
                    }
                )

            # sort so genes end up being in correct order for overwriting
            correctedGenes <-
                correctedGenes[order(
                    correctedGenes$isoform_id,
                    correctedGenes$gene_id,
                    correctedGenes$sample_1,
                    correctedGenes$sample_2
                ), -1] # -1 removes the index created by ldply
            # overwrite problematic genes
            isoformData[geneNameIndexesData, c(
                'gene_id',
                'sample_1',
                'sample_2',
                'gene_value_1',
                'gene_value_2',
                'gene_stderr_1',
                'gene_stderr_2',
                'gene_log2_fold_change',
                'gene_p_value',
                'gene_q_value',
                'gene_significant',
                'iso_status',
                'iso_value_1',
                'iso_value_2'
            )] <- correctedGenes[, -1] # -1 removes the isoform id

            if (!quiet) {
                message(
                    paste(
                        "Cufflinks annotation problem was fixed for",
                        length(geneNameSplit),
                        "Cuff_genes",
                        sep = ' '
                    )
                )
            }
        } else {
            if (!quiet) {
                message(
                    "No instances of a Cufflinks annotation problem found - no changes were made"
                )
            }
        }
    }

    ### Calculate IF values
    isoformData$IF1 <-
        isoformData$iso_value_1 / isoformData$gene_value_1
    isoformData$IF2 <-
        isoformData$iso_value_2 / isoformData$gene_value_2
    isoformData$dIF <- isoformData$IF2 - isoformData$IF1


    ### Add q-values
    if (addCufflinksSwichTest) {
        if (!quiet) {
            message("Extracting analysis of alternative splicing...\n")
        }

        ### Extract cufflinks switch test
        myCuffGeneSet <- suppressMessages(
            cummeRbund::getGenes(cuffDB, featureNames(genes(cuffDB)))
        )
        cuffSplicing <- cummeRbund::diffData(
            cummeRbund::splicing(myCuffGeneSet)
        )

        isoformData$isoform_switch_q_value <- NA

        if (nrow(cuffSplicing)) {
            isoformData$gene_switch_q_value <- cuffSplicing$q_value[match(
                paste(
                    isoformData$gene_id,
                    isoformData$sample_1,
                    isoformData$sample_2,
                    sep = '_'
                ),
                paste(
                    cuffSplicing$gene_id,
                    cuffSplicing$sample_1,
                    cuffSplicing$sample_2,
                    sep = '_'
                )
            )]
        } else {
            isoformData$gene_switch_q_value <- NA
        }
    } else {
        ### Add collumns for isoform switch analysis results
        isoformData$isoform_switch_q_value <- NA
        isoformData$gene_switch_q_value <- NA
    }

    if (!quiet) {
        message("Making IsoformSwitchAanalyzeRlist...\n")
    }

    ### Reorder a bit
    isoformData <- isoformData[, c(1, 4, 2:3, 5:ncol(isoformData))]
    colnames(isoformData)[3:4] <- c('condition_1', 'condition_2')

    # Create GRanges for exon features
    exonFeatures <- GRanges(
        seqnames = isoformFeatures$"seqnames",
        strand = isoformFeatures$"strand",
        ranges = IRanges(start = isoformFeatures$"start",
                         end = isoformFeatures$"end"),
        isoformFeatures[, c("isoform_id", "gene_id")]
    )

    ### Extract run info
    # cufflinks version
    cuffVersion <- cummeRbund::runInfo(cuffDB)$value[2]

    ### Check cufflinks version
    checkVersionFail <- function(versionVector, minVersionVector) {
        for (i in seq_along(versionVector)) {
            if (versionVector[i] > minVersionVector[i]) {
                return(FALSE)
            }
            if (versionVector[i] < minVersionVector[i]) {
                return(TRUE)
            }
        }
        return(FALSE)
    }

    cuffVersionDeconstructed <-
        as.integer(strsplit(cuffVersion, '\\.')[[1]])
    if (checkVersionFail(cuffVersionDeconstructed, c(2, 2, 1))) {
        warning(
            paste(
                'The version of cufflinks/cuffdiff you have used is outdated',
                '. An error in the estimations of standard deviations was not',
                'corrected untill cufflinks 2.2.1. Since this detection of',
                'isoform switches using this R package relies',
                'on this standard error estimat',
                'we do not premit switch detection',
                '(using the detectIsoformSwitches() ) with the data generated here.',
                'If you want to use this R pacakge please',
                'upgrade cufflinks/cuffdiff to version >=2.2.1 or newer and try again.',
                sep = ' '
            )
        )
    }

    # replicate numbers
    repInfo <- cummeRbund::replicates(cuffDB)
    nrRep <- table(repInfo$sample_name)
    nrRep <-
        data.frame(
            condition = names(nrRep),
            nrReplicates = as.vector(nrRep),
            row.names = NULL,
            stringsAsFactors = FALSE
        )

    ### Design matrix
    designMatrix <- repInfo[, c('rep_name', 'sample_name')]
    colnames(designMatrix) <- c('sampleID', 'condition')

    # Return SpliceRList
    switchAnalyzeRlist <- createSwitchAnalyzeRlist(
        isoformFeatures = isoformData,
        exons = exonFeatures,
        designMatrix = designMatrix,
        isoformCountMatrix = isoRepExp,
        sourceId = paste("cufflinks", cuffVersion , sep = '_')
    )

    if (addCufflinksSwichTest & nrow(cuffSplicing)) {
        switchAnalyzeRlist$isoformSwitchAnalysis <- cuffSplicing
    }

    if (!quiet) {
        message('Done')
    }
    return(switchAnalyzeRlist)
}

importCufflinksFiles <- function(
    pathToGTF,
    pathToGeneDEanalysis,
    pathToIsoformDEanalysis,
    pathToGeneFPKMtracking,
    pathToIsoformFPKMtracking,
    pathToIsoformReadGroupTracking,
    pathToSplicingAnalysis = NULL,
    pathToReadGroups,
    pathToRunInfo,
    fixCufflinksAnnotationProblem = TRUE,
    quiet = FALSE
) {
    ### Test that files exist
    if (TRUE) {
        # pathToGTF
        if (!file.exists(pathToGTF)) {
            stop('The \'pathToGTF\' argument does not point to an acutal file')
        }
        # DE
        if (!file.exists(pathToGeneDEanalysis))    {
            stop('The \'pathToGeneDEanalysis\' argument does not point to an acutal file')
        }
        if (!file.exists(pathToIsoformDEanalysis)) {
            stop('The \'pathToIsoformDEanalysis\' argument does not point to an acutal file')
        }
        # Tracking
        if (!file.exists(pathToGeneFPKMtracking))    {
            stop('The \'pathToGeneFPKMtracking\' argument does not point to an acutal file')
        }
        if (!file.exists(pathToIsoformFPKMtracking)) {
            stop(
                'The \'pathToIsoformFPKMtracking\' argument does not point to an acutal file'
            )
        }
        if (!file.exists(pathToIsoformReadGroupTracking)) {
            stop(
                'The \'pathToIsoformReadGroupTracking\' argument does not point to an acutal file'
            )
        }
        # splicing
        if (!is.null(pathToSplicingAnalysis)) {
            if (!file.exists(pathToSplicingAnalysis)) {
                stop(
                    'The \'pathToSplicingAnalysis\' argument does not point to an acutal file'
                )
            }
        }
        # info
        if (!file.exists(pathToReadGroups)) {
            stop('The \'pathToReadGroups\' argument does not point to an acutal file')
        }
        if (!file.exists(pathToRunInfo))    {
            stop('The \'pathToRunInfo\' argument does not point to an acutal file')
        }
    }

    ### Import the supplied files (not gtf)
    if (TRUE) {
        if (!quiet) {
            message("Loading genes and isoforms...")
        }
        geneDiffanalysis     <-
            read.table(
                file = pathToGeneDEanalysis,
                header = TRUE,
                stringsAsFactors = FALSE,
                sep = '\t'
            )
        isoformDiffanalysis  <-
            read.table(
                file = pathToIsoformDEanalysis,
                header = TRUE,
                stringsAsFactors = FALSE,
                sep = '\t'
            )

        geneAnnotation       <-
            read.table(
                file = pathToGeneFPKMtracking,
                header = TRUE,
                stringsAsFactors = FALSE,
                sep = '\t'
            )
        isoformAnnotation    <-
            read.table(
                file = pathToIsoformFPKMtracking,
                header = TRUE,
                stringsAsFactors = FALSE,
                sep = '\t'
            )
        isoRepExp        <-
            read.table(
                file = pathToIsoformReadGroupTracking,
                header = TRUE,
                stringsAsFactors = FALSE,
                sep = '\t'
            )

        cuffSplicing         <-
            read.table(
                file = pathToSplicingAnalysis,
                header = TRUE,
                stringsAsFactors = FALSE,
                sep = '\t'
            )

        readGroup <-
            read.table(
                file = pathToReadGroups,
                header = TRUE,
                stringsAsFactors = FALSE,
                sep = '\t'
            )
        runInfo   <-
            read.table(
                file = pathToRunInfo,
                header = TRUE,
                stringsAsFactors = FALSE,
                sep = '\t'
            )
    }

    ### "Test" that the data.files are what they are supposed to be
    if (TRUE) {
        ### gene diff analysis
        q1 <-
            !all(
                colnames(geneDiffanalysis) %in% c(
                    "test_id",
                    "gene_id",
                    "gene",
                    "locus",
                    "sample_1",
                    "sample_2",
                    "status",
                    "value_1",
                    "value_2",
                    "log2.fold_change.",
                    "test_stat",
                    "p_value",
                    "q_value",
                    "significant"
                )
            )
        q2 <-
            sum(grepl(
                'XLOC', geneDiffanalysis$test_id
            )) != nrow(geneDiffanalysis)
        if (q1 | q2) {
            stop(paste(
                'The file supplied to pathToGeneDEanalysis does not appear',
                'to be the result of the CuffDiff gene expression analysis.'
            ))
        }
        ### transcript diff analysis
        q1 <-
            !all(
                colnames(isoformDiffanalysis) %in% c(
                    "test_id",
                    "gene_id",
                    "gene",
                    "locus",
                    "sample_1",
                    "sample_2",
                    "status",
                    "value_1",
                    "value_2",
                    "log2.fold_change.",
                    "test_stat",
                    "p_value",
                    "q_value",
                    "significant"
                )
            )
        q2 <-
            sum(grepl(
                'TCONS', isoformDiffanalysis$test_id
            )) != nrow(isoformDiffanalysis)
        if (q1 | q2) {
            stop(paste(
                'The file supplied to isoformDiffanalysis does not appear to',
                'be the result of the CuffDiff transcript expression analysis.'
            ))
        }

        ### gene annoation
        q1 <-
            !all(
                colnames(geneAnnotation)[1:8] %in% c(
                    "tracking_id",
                    "class_code",
                    "nearest_ref_id",
                    "gene_id",
                    "gene_short_name",
                    "tss_id",
                    "locus",
                    "length"
                )
            )
        q2 <-
            sum(grepl(
                'XLOC', geneAnnotation$tracking_id
            )) != nrow(geneAnnotation)
        if (q1 | q2) {
            stop(paste(
                'The file supplied to geneAnnotation does not appear to be the',
                'gene FPKM traccking of the CuffDiff gene FPKM trascking analysis.'
            ))
        }
        ### transcript annoation
        q1 <-
            !all(
                colnames(isoformAnnotation)[1:8] %in% c(
                    "tracking_id",
                    "class_code",
                    "nearest_ref_id",
                    "gene_id",
                    "gene_short_name",
                    "tss_id",
                    "locus",
                    "length"
                )
            )
        q2 <-
            sum(grepl(
                'TCONS', isoformAnnotation$tracking_id
            )) != nrow(isoformAnnotation)
        if (q1 | q2) {
            stop(paste(
                'The file supplied to isoformAnnotation does not appear to be',
                'the isoform FPKM tracking of the CuffDiff transcript analysis.'
            ))
        }

        ### rep expression
        q1 <-
            !all(
                colnames(isoRepExp)[1:4] %in% c(
                    "tracking_id", "condition", "replicate", "raw_frags"
                )
            )
        q2 <-
            sum(grepl('TCONS', isoRepExp$tracking_id)) != nrow(isoRepExp)
        if (q1 | q2) {
            stop(paste(
                'The file supplied to pathToIsoformCountTracking does not',
                'appear to be the isoform count tracking of the CuffDiff',
                'transcript analysis.'
            ))
        }

        ### splicing analysis
        q1 <-
            !all(
                colnames(cuffSplicing) %in% c(
                    "test_id",
                    "gene_id",
                    "gene",
                    "locus",
                    "sample_1",
                    "sample_2",
                    "status",
                    "value_1",
                    "value_2",
                    "sqrt.JS.",
                    "test_stat",
                    "p_value",
                    "q_value",
                    "significant"
                )
            )
        q2 <-
            sum(grepl('TSS', cuffSplicing$test_id)) != nrow(cuffSplicing)
        if (q1 | q2) {
            stop(
                'The file supplied to cuffSplicing does not appear to be the',
                'result of the CuffDiff differential analysis of alternative splicing.'
            )
        }

        ### Read grous
        q1 <-
            !all(
                colnames(readGroup) %in% c(
                    "file",
                    "condition",
                    "replicate_num",
                    "total_mass",
                    "norm_mass",
                    "internal_scale",
                    "external_scale"
                )
            )
        q2 <-
            !all(readGroup$condition %in% unique(
                unlist(geneDiffanalysis[, c('sample_1', 'sample_2')]))
            )
        if (q1 | q2) {
            stop(paste(
                'The file supplied to readGroup does not appear to be the',
                'pathToReadGroups of the CuffDiff transcript analysis.'
            ))
        }

        ### Run info
        q1 <- !all(colnames(runInfo) %in% c("param", 'value'))
        q2 <-
            !all(
                runInfo$param %in% c(
                    "cmd_line",
                    "version",
                    "SVN_revision",
                    "boost_version"
                )
            )
        if (q1 | q2) {
            stop(paste(
                'The file supplied to runInfo does not appear to be',
                'the runInfo of the CuffDiff transcript analysis.'
            ))
        }

    }

    ### Massage and merge gene and isoform annoation and DE analysis
    if (TRUE) {
        if (!quiet) {
            message("Merging gene and isoform expression...")
        }
        ### Massage data frames
        if (TRUE) {
            # gene
            geneDiffanalysis <-
                geneDiffanalysis[, -which(
                    colnames(geneDiffanalysis) %in%
                        c('test_id', 'gene', 'locus', 'test_stat')
                )]
            colnames(geneDiffanalysis)[4:ncol(geneDiffanalysis)] <-
                paste(
                    "gene_",
                    colnames(geneDiffanalysis)[4:ncol(geneDiffanalysis)] ,
                    sep = "") # add gene to the colnames so they can be destinquished from the gene diff data
            colnames(geneDiffanalysis) <- gsub(
                    'gene_log2.fold_change.',
                    'gene_log2_fold_change',
                    colnames(geneDiffanalysis)
                )

            # info
            colnames(isoformAnnotation)[1] <- 'isoform_id'
            isoformAnnotation2 <-
                isoformAnnotation[, na.omit(match(
                    c(
                        'isoform_id',
                        'gene_id',
                        'gene_short_name',
                        'nearest_ref_id',
                        'class_code',
                        'tss_id',
                        'CDS_id',
                        'length',
                        'locus'
                    ),
                    colnames(isoformAnnotation)
                ))]

            colnames(isoformAnnotation2)[which(
                colnames(isoformAnnotation2) == 'gene_short_name'
            )] <- 'gene_name'

            # iso
            isoformDiffanalysis <- isoformDiffanalysis[, which(
                !colnames(isoformDiffanalysis) %in%
                    c('gene_id', 'gene', 'test_stat')
            )]
            colnames(isoformDiffanalysis)[5:ncol(isoformDiffanalysis)] 	<-
                paste(
                    "iso_",
                    colnames(isoformDiffanalysis)[5:ncol(isoformDiffanalysis)],
                    sep = "") # add gene to the colnames so they can be destinquished from the gene diff data
            colnames(isoformDiffanalysis)[1] <- 'isoform_id'
            colnames(isoformDiffanalysis) <-
                gsub(
                    'iso_log2.fold_change.',
                    'iso_log2_fold_change',
                    colnames(isoformDiffanalysis)
                )

            # rep expression
            isoRepExp <-
                isoRepExp[, c("tracking_id",
                              "condition",
                              "replicate",
                              "raw_frags")]
            colnames(isoRepExp)[1] <- 'isoform_id'
            isoRepExp$rep <-
                paste(isoRepExp$condition, isoRepExp$replicate, sep = '_')
            isoRepExp <-
                reshape2::dcast(data = isoRepExp,
                                isoform_id ~ rep,
                                value.var = 'raw_frags')

        }

        ### Extract standard error
        if (TRUE) {
            ### Tjek if the Isoform CI collums are switches
            ciLowColumn <-
                which(grepl('_conf_lo', colnames(isoformAnnotation)))[1]
            ciHighColumn <-
                which(grepl('_conf_hi', colnames(isoformAnnotation)))[1]

            if (all(isoformAnnotation[, ciHighColumn] >= isoformAnnotation[, ciLowColumn])) {
                highString <- '_conf_hi'
            } else {
                highString <- '_conf_lo'
            }

            ### extract isoform sddev from CI
            # fpkm
            isoformFPKM <-
                isoformAnnotation[, which(grepl(
                    'isoform_id|_FPKM', colnames(isoformAnnotation)
                ))]
            isoformFPKM <- melt(isoformFPKM, id.vars = 'isoform_id')
            isoformFPKM$variable <-
                gsub('_FPKM$', '', isoformFPKM$variable)
            colnames(isoformFPKM)[3] <- 'expression'
            # ci high
            isoformFPKMciHi <-
                isoformAnnotation[, which(grepl(
                    paste('isoform_id|', highString, sep = ''),
                    colnames(isoformAnnotation)
                ))]
            isoformFPKMciHi <-
                melt(isoformFPKMciHi, id.vars = 'isoform_id')
            isoformFPKMciHi$variable <-
                gsub(highString, '', isoformFPKMciHi$variable)
            colnames(isoformFPKMciHi)[3] <- 'ci_hi'
            # stderr
            isoformFPKMcombined <-
                merge(isoformFPKM,
                      isoformFPKMciHi,
                      by = c('isoform_id', 'variable'))
            isoformFPKMcombined$iso_stderr <-
                (isoformFPKMcombined$ci_hi - isoformFPKMcombined$expression) / 2 # How it's done in cufflinks source code
            isoformFPKMcombined$expression <- NULL
            isoformFPKMcombined$ci_hi <- NULL
            colnames(isoformFPKMcombined) <-
                c('isoform_id', 'sample_name', 'iso_stderr')

            ### Tjek if the gene CI collums are switches
            ciLowColumn <-
                which(grepl('_conf_lo', colnames(geneAnnotation)))[1]
            ciHighColumn <-
                which(grepl('_conf_hi', colnames(geneAnnotation)))[1]

            if (all(
                geneAnnotation[, ciHighColumn] >= geneAnnotation[, ciLowColumn])
            ) {
                highString <- '_conf_hi'
            } else {
                highString <- '_conf_lo'
            }

            ### extract gene sddev from CI
            # fpkm
            geneFPKM <- geneAnnotation[, which(grepl(
                    'tracking_id|_FPKM', colnames(geneAnnotation)
                ))]
            geneFPKM <- melt(geneFPKM, id.vars = 'tracking_id')
            geneFPKM$variable <-
                gsub('_FPKM$', '', geneFPKM$variable)
            colnames(geneFPKM)[3] <- 'expression'
            # ci high
            geneFPKMciHi <- geneAnnotation[, which(grepl(
                    paste('tracking_id|', highString, sep = ''),
                    colnames(geneAnnotation)
                ))]
            geneFPKMciHi <- melt(geneFPKMciHi, id.vars = 'tracking_id')
            geneFPKMciHi$variable <- gsub(highString, '', geneFPKMciHi$variable)
            colnames(geneFPKMciHi)[3] <- 'ci_hi'
            # stderr
            geneFPKMcombined <-
                merge(geneFPKM,
                      geneFPKMciHi,
                      by = c('tracking_id', 'variable'))
            geneFPKMcombined$iso_stderr <-
                (geneFPKMcombined$ci_hi - geneFPKMcombined$expression) / 2 # how it's done in cufflinks sourece code
            geneFPKMcombined$expression <- NULL
            geneFPKMcombined$ci_hi <- NULL
            colnames(geneFPKMcombined) <-
                c('gene_id', 'sample_name', 'gene_stderr')


            ## Merge stderr with DE analysis
            isoformDiffanalysis <-
                merge(
                    isoformDiffanalysis,
                    isoformFPKMcombined,
                    by.x = c('isoform_id', 'sample_2'),
                    by.y = c('isoform_id', 'sample_name')
                )
            colnames(isoformDiffanalysis)[which( grepl(
                'iso_stderr', colnames(isoformDiffanalysis))
            )] <- 'iso_stderr_2'

            isoformDiffanalysis <- merge(
                    isoformDiffanalysis,
                    isoformFPKMcombined,
                    by.x = c('isoform_id', 'sample_1'),
                    by.y = c('isoform_id', 'sample_name')
                )
            colnames(isoformDiffanalysis)[which(grepl(
                'iso_stderr$',
                colnames(isoformDiffanalysis),
                perl = TRUE
            ))] <- c('iso_stderr_1')

            isoformDiffanalysis <-
                isoformDiffanalysis[, c(
                    'isoform_id',
                    'sample_1',
                    'sample_2',
                    'iso_status',
                    'iso_value_1',
                    'iso_value_2',
                    'iso_stderr_1',
                    'iso_stderr_2',
                    'iso_log2_fold_change',
                    'iso_p_value',
                    'iso_q_value',
                    'iso_significant'
                )]

            ### Extract and add gene stderr
            geneDiffanalysis <-
                merge(
                    geneDiffanalysis,
                    geneFPKMcombined,
                    by.x = c('gene_id', 'sample_2'),
                    by.y = c('gene_id', 'sample_name')
                )
            colnames(geneDiffanalysis)[ which(grepl(
                'gene_stderr', colnames(geneDiffanalysis)
            ))] <- 'gene_stderr_2'

            geneDiffanalysis <- merge(
                    geneDiffanalysis,
                    geneFPKMcombined,
                    by.x = c('gene_id', 'sample_1'),
                    by.y = c('gene_id', 'sample_name')
                )
            colnames(geneDiffanalysis)[which(grepl(
                'gene_stderr$',
                colnames(geneDiffanalysis),
                perl = TRUE
            ))] <- c('gene_stderr_1')

            geneDiffanalysis <-
                geneDiffanalysis[, c(
                    'gene_id',
                    'sample_1',
                    'sample_2',
                    'gene_status',
                    'gene_value_1',
                    'gene_value_2',
                    'gene_stderr_1',
                    'gene_stderr_2',
                    'gene_log2_fold_change',
                    'gene_p_value',
                    'gene_q_value',
                    'gene_significant'
                )]

        }

        ### Merge data
        if (TRUE) {
            ### Meger gene DE and annotation data
            isoformData <-
                merge(isoformAnnotation2, geneDiffanalysis, by = 'gene_id')

            ### Merge with iso DE
            isoformData <-
                merge(
                    isoformData,
                    isoformDiffanalysis,
                    by = c('isoform_id', 'sample_1', 'sample_2')
                )

            ### Massage again
            colnames(isoformData)[which(
                colnames(isoformData) == 'tss_id'
            )] <- 'TSS_group_id'

        }

    }


    ### Obtain transcript structure information
    if (TRUE) {
        if (!quiet) {
            message(paste(
                'Importing transcript structure. This can take a',
                'while (since a GTF is a unstructured file)...'
            ))
        }
        ### Import file
        exonFeatures <-  rtracklayer::import(pathToGTF, format = 'gtf')
        if (length(exonFeatures) == 0)
            stop("No exon information extracted from GTF")

        ### Filter for what is needed
        exonFeatures <- exonFeatures[
            which(tolower(exonFeatures$type) == 'exon'),
            c('gene_id', 'transcript_id')
            ]

        ### rename
        colnames(exonFeatures@elementMetadata) <- gsub(
            'transcript_id', 'isoform_id',
            colnames(exonFeatures@elementMetadata)
        )
    }

    ### Check it is the same transcripts in transcript structure and expression info
    if (TRUE) {
        myUnion     <-
            unique(c(
                isoformData$isoform_id,
                exonFeatures$isoform_id,
                isoRepExp$isoform_id
            ))
        myIntersect <- intersect(
            intersect(isoformData$isoform_id, exonFeatures$isoform_id),
            isoRepExp$isoform_id
        )

        # If there are descripencies
        if (length(myUnion) != length(myIntersect)) {
            isoformData <- isoformData[which(
                isoformData$isoform_id     %in% myIntersect), ]
            exonFeatures <- exonFeatures[which(
                    exonFeatures$isoform_id    %in% myIntersect), ]
            isoRepExp <- isoRepExp[which(
                isoRepExp$isoform_id       %in% myIntersect), ]

            if (!quiet) {
                message(
                    paste(
                        'There were discrepencies between the GTF and the',
                        'expression analysis files. To solve this',
                        abs(length(myUnion) - length(myIntersect)) ,
                        'transcripts were removed.',
                        sep = ' '
                    )
                )
            }
        }
    }

    ### Fix to correct for Cufflinks annotation problem where cufflinks assignes
    # transcripts from several annotated genes to 1 cuffgene
    if (fixCufflinksAnnotationProblem) {
        if (!quiet) {
            message("Analyzing cufflinks annotation problem...")
        }

        geneName <- unique(isoformData[, c('gene_id', 'gene_name')])
        geneNameSplit <-
            split(geneName$gene_name , f = geneName$gene_id)
        # remove all unique
        geneNameSplit <-
            geneNameSplit[which(sapply(geneNameSplit, function(x)
                length(unique(x))) > 1)]

        if (length(geneNameSplit) > 0) {
            # if there are any problems
            if (!quiet) {
                message("Fixing cufflinks annotation problem...")
            }
            #get indexes of those affected
            geneNameIndexesData     <-
                which(isoformData$gene_id %in% names(geneNameSplit))
            geneNameIndexesFeatures <-
                which(exonFeatures$spliceR.gene_id %in% names(geneNameSplit))

            # combine names of cuffgenes and gene short name
            isoformData$gene_id[geneNameIndexesData]              <-
                paste(isoformData$gene_id[geneNameIndexesData]        ,
                      isoformData$gene_name[geneNameIndexesData],
                      sep = ':')
            exonFeatures$spliceR.gene_id[geneNameIndexesFeatures] <-
                paste(exonFeatures$spliceR.gene_id[geneNameIndexesFeatures],
                      exonFeatures$spliceR.gene_name[geneNameIndexesFeatures],
                      sep = ':')

            ## Correct gene expression levels and differntial analysis
            problematicGenes <-
                isoformData[geneNameIndexesData, c(
                    'isoform_id',
                    'gene_id',
                    'sample_1',
                    'sample_2',
                    'gene_value_1',
                    'gene_value_2',
                    'gene_stderr_1',
                    'gene_stderr_2',
                    'gene_log2_fold_change',
                    'gene_p_value',
                    'gene_q_value',
                    'gene_significant',
                    'iso_status',
                    'iso_value_1',
                    'iso_value_2'
                )]
            problematicGenesSplit <-
                split(problematicGenes, f = problematicGenes[
                    ,c('gene_id', 'sample_1', 'sample_2')], drop =TRUE)

            correctedGenes <-
                ldply(
                    problematicGenesSplit,
                    .fun = function(df) {
                        # df <- problematicGenesSplit[[1]]
                        df$gene_value_1 <- sum(df$iso_value_1)
                        df$gene_value_2 <- sum(df$iso_value_2)
                        df$gene_stderr_1 <- NA
                        df$gene_stderr_2 <- NA
                        df$gene_log2_fold_change <- log2(
                            (df$gene_value_2[2] + 0.0001) /
                                (df$gene_value_1[1] + 0.0001)
                        )
                        df$gene_p_value <- 1
                        df$gene_q_value <- 1
                        df$gene_significant <- 'no'
                        df$iso_status <- 'NOTEST'
                        return(df)
                    }
                )

            # sort so genes end up being in correct order for overwriting
            correctedGenes <-
                correctedGenes[order(
                    correctedGenes$isoform_id,
                    correctedGenes$gene_id,
                    correctedGenes$sample_1,
                    correctedGenes$sample_2
                ), -1] # -1 removes the index created by ldply
            # overwrite problematic genes
            isoformData[geneNameIndexesData, c(
                'gene_id',
                'sample_1',
                'sample_2',
                'gene_value_1',
                'gene_value_2',
                'gene_stderr_1',
                'gene_stderr_2',
                'gene_log2_fold_change',
                'gene_p_value',
                'gene_q_value',
                'gene_significant',
                'iso_status',
                'iso_value_1',
                'iso_value_2'
            )] <- correctedGenes[, -1] # -1 removes the isoform id

            if (!quiet) {
                message(
                    paste(
                        "Cufflinks annotation problem was fixed for",
                        length(geneNameSplit),
                        "Cuff_genes",
                        sep = ' '
                    )
                )
            }
        } else {
            if (!quiet) {
                message(paste(
                    'No instances of a Cufflinks annotation',
                    'problem found - no changes were made'
                ))
            }
        }
    } # end of fix cufflinks annotatopn problem

    if (!quiet) {
        message("Making IsoformSwitchAanalyzeRlist...\n")
    }

    ### Calculate IF values
    isoformData$IF1 <-
        isoformData$iso_value_1 / isoformData$gene_value_1
    isoformData$IF2 <-
        isoformData$iso_value_2 / isoformData$gene_value_2
    isoformData$dIF <- isoformData$IF2 - isoformData$IF1

    ### Add q-values
    if (!is.null(pathToSplicingAnalysis)) {
        ### Add cufflinks analysis isoform switch analysis results
        isoformData$isoform_switch_q_value <- NA
        isoformData$gene_switch_q_value <-
            cuffSplicing$q_value[match(
                paste(
                    isoformData$gene_id,
                    isoformData$sample_1,
                    isoformData$sample_2,
                    sep = '_'
                ),
                paste(
                    cuffSplicing$gene_id,
                    cuffSplicing$sample_1,
                    cuffSplicing$sample_2,
                    sep = '_'
                )
            )]
    } else {
        ### Add collumns for isoform switch analysis results
        isoformData$isoform_switch_q_value <- NA
        isoformData$gene_switch_q_value <- NA
    }

    ### Reorder a bit
    isoformData <- isoformData[, c(1, 4, 2:3, 5:ncol(isoformData))]
    colnames(isoformData)[3:4] <- c('condition_1', 'condition_2')

    ### Extract run info
    # cufflinks version
    cuffVersion <- runInfo$value[2]

    ### Check cufflinks version
    checkVersionFail <- function(versionVector, minVersionVector) {
        for (i in seq_along(versionVector)) {
            if (versionVector[i] > minVersionVector[i]) {
                return(FALSE)
            }
            if (versionVector[i] < minVersionVector[i]) {
                return(TRUE)
            }
        }
        return(FALSE)
    }

    cuffVersionDeconstructed <-
        as.integer(strsplit(cuffVersion, '\\.')[[1]])
    if (checkVersionFail(cuffVersionDeconstructed, c(2, 2, 1))) {
        warning(
            paste(
                'The version of cufflinks/cuffdiff you have used is outdated',
                '. An error in the estimations of standard deviations was not',
                'corrected untill cufflinks 2.2.1. Since this detection of',
                'isoform switches using this R package relies',
                'on this standard error estimat',
                'we do not premit switch detection',
                '(using the detectIsoformSwitches() ) with the data generated here.',
                'If you want to use this R pacakge please',
                'upgrade cufflinks/cuffdiff to version >=2.2.1 or newer and try again.',
                sep = ' '
            )
        )
    }

    # replicate numbers
    nrRep <- table(readGroup$condition)
    nrRep <-
        data.frame(
            condition = names(nrRep),
            nrReplicates = as.vector(nrRep),
            row.names = NULL,
            stringsAsFactors = FALSE
        )

    ### Design matrix
    readGroup$sample_name <-
        paste0(readGroup$condition, '_', readGroup$replicate_num)
    designMatrix <- readGroup[, c('sample_name', 'condition')]
    colnames(designMatrix) <- c('sampleID', 'condition')

    # Return SpliceRList
    switchAnalyzeRlist <- createSwitchAnalyzeRlist(
        isoformFeatures = isoformData,
        exons = exonFeatures,
        designMatrix = designMatrix,
        isoformCountMatrix = isoRepExp,
        sourceId = paste("cufflinks", cuffVersion , sep = '_')
    )

    if (!is.null(pathToSplicingAnalysis) & nrow(cuffSplicing)) {
        switchAnalyzeRlist$isoformSwitchAnalysis <- cuffSplicing
    }

    return(switchAnalyzeRlist)
}

importGTF <- function(
    pathToGTF,
    addAnnotatedORFs = FALSE,
    onlyConsiderFullORF = FALSE,
    removeNonConvensionalChr = FALSE,
    PTCDistance = 50,
    quiet = FALSE
) {
    # Read in from GTF file and create Rdata file for easy loading
    if (!quiet) {
        message('importing GTF (this may take a while)')
    }
    mfGTF <- rtracklayer::import(pathToGTF)

    ### tjeck GTF
    if (!all(c('transcript_id', 'gene_id') %in% colnames(mfGTF@elementMetadata))) {
        collumnsMissing <- paste(
            c('transcript_id', 'gene_id')[which(
                !c('transcript_id', 'gene_id') %in%
                    colnames(mfGTF@elementMetadata)
            )], collapse = ', ')
        stop(
            paste(
                'The GTF file must contain the folliwing collumns',
                '\'transcript_id\', \'gene_id\' and.',
                collumnsMissing,
                'is missing.',
                sep = ' '
            )
        )
    }

    ### Reduce if nessesary
    if (removeNonConvensionalChr) {
        mfGTF <- mfGTF[which(!grepl('_', as.character(mfGTF@seqnames))), ]

        if (length(mfGTF) == 0) {
            stop('No exons were left after filtering',
                 'with \'removeNonConvensionalChr\'.')
        }

        seqlevels(mfGTF) <- as.character(mfGTF@seqnames@values)
    }

    ### Make annoation
    if (!quiet) {
        message('converting GTF to switchAnalyzeRlist')
    }
    exonAnoationIndex <- which(mfGTF$type == 'exon')

    colsToExtract <- c('transcript_id', 'gene_id', 'gene_name')
    myIso <-
        as.data.frame(unique(mfGTF@elementMetadata[
            exonAnoationIndex,
            na.omit(match(colsToExtract, colnames(mfGTF@elementMetadata)))]
        ))

    if (is.null(myIso$gene_name)) {
        myIso$gene_name <- NA
    }

    myIsoAnot <- data.frame(
        isoform_id = myIso$transcript_id,
        gene_id = myIso$gene_id,
        condition_1 = "plaseholder1",
        condition_2 = "plaseholder2",
        gene_name = myIso$gene_name,
        class_code = NA,
        gene_value_1 = 0,
        gene_value_2 = 0,
        gene_stderr_1 = NA,
        gene_stderr_2 = NA,
        gene_log2_fold_change = 0,
        gene_p_value = 1,
        gene_q_value = 1,
        iso_value_1 = 0,
        iso_value_2 = 0,
        iso_stderr_1 = NA,
        iso_stderr_2 = NA,
        iso_log2_fold_change = 0,
        iso_p_value = 1,
        iso_q_value = 1,
        IF1 = NA,
        IF2 = NA,
        dIF = NA,
        isoform_switch_q_value = NA,
        gene_switch_q_value = NA,
        stringsAsFactors = FALSE
    )

    ### Add CDS annoation from GTF file inc convertion to transcript coordinats
    if (addAnnotatedORFs) {
        # test whether any CDS are found
        if (any(mfGTF$type == 'CDS')) {
            message('converting annotated CDSs')
            myCDS <-
                sort(mfGTF[which(mfGTF$type == 'CDS'), 'transcript_id'])
            myCDSedges <-
                suppressMessages(unlist(range(
                    split(myCDS[, 0], f = myCDS$transcript_id)
                )))  # Extract EDGEs
            myCDSedges$id <- names(myCDSedges)
            names(myCDSedges) <- NULL

            if (onlyConsiderFullORF) {
                fullyAnnoated <-
                    as.data.frame(sort(
                        mfGTF[which(
                            mfGTF$type %in% c('start_codon', 'stop_codon')),
                            c('transcript_id', 'type')]))
                fullyAnnoatedSplit <-
                    split(as.character(fullyAnnoated$type),
                          f = fullyAnnoated$transcript_id)
                fullyAnnoatedCount <-
                    sapply(fullyAnnoatedSplit, function(x)
                        length(unique(x)))
                toKeep <-
                    names(fullyAnnoatedCount[which(fullyAnnoatedCount == 2)])


                myCDSedges <-
                    myCDSedges[which(myCDSedges$id %in% toKeep), ]
            }

            ### Extract Exons
            localExons <- mfGTF[exonAnoationIndex, 'transcript_id']
            localExons <-
                localExons[which(
                    as.character(localExons@strand) %in% c('+', '-')), ]
            localExons <-
                localExons[which(localExons$transcript_id %in% myCDSedges$id), ]

            localExons <-
                localExons[order(localExons$transcript_id,
                                 start(localExons),
                                 end(localExons)), ]
            localExons$exon_id <-
                paste('exon_', 1:length(localExons), sep = '')

            ### Extract strand specific ORF info
            cds <- as.data.frame(myCDSedges)
            # start
            plusIndex <- which(cds$strand == '+')
            annoatedStartGRangesPlus <-
                GRanges(
                    cds$seqnames[plusIndex],
                    IRanges(
                        start = cds$start[plusIndex],
                        end = cds$start[plusIndex]),
                    strand = cds$strand[plusIndex],
                    id = cds$id[plusIndex]
                )
            minusIndex <- which(cds$strand == '-')
            annoatedStartGRangesMinus <-
                GRanges(
                    cds$seqnames[minusIndex],
                    IRanges(
                        start = cds$end[minusIndex],
                        end = cds$end[minusIndex]),
                    strand = cds$strand[minusIndex],
                    id = cds$id[minusIndex]
                )

            annoatedStartGRanges <-
                c(annoatedStartGRangesPlus,
                  annoatedStartGRangesMinus)
            annoatedStartGRanges$orf_id <-
                paste('cds_', 1:length(annoatedStartGRanges), sep = '')

            # end
            annoatedEndGRangesPlus  <-
                GRanges(
                    cds$seqnames[plusIndex],
                    IRanges(
                        start = cds$end[plusIndex],
                        end = cds$end[plusIndex]),
                    strand = cds$strand[plusIndex],
                    id = cds$id[plusIndex]
                )
            annoatedEndGRangesMinus <-
                GRanges(
                    cds$seqnames[minusIndex],
                    IRanges(
                        start = cds$start[minusIndex],
                        end = cds$start[minusIndex]),
                    strand = cds$strand[minusIndex],
                    id = cds$id[minusIndex]
                )

            annoatedEndGRanges <-
                c(annoatedEndGRangesPlus, annoatedEndGRangesMinus)
            annoatedEndGRanges$orf_id <-
                paste('stop_', 1:length(annoatedEndGRanges), sep = '')

            # combine
            annotatedORFGR <-
                c(annoatedStartGRanges, annoatedEndGRanges)


            ### Idenetify overlapping CDS and exons as well as the annoate transcript id
            suppressWarnings(overlappingAnnotStart <-
                                 as.data.frame(
                                     findOverlaps(
                                         query = localExons,
                                         subject = annotatedORFGR,
                                         ignore.strand = FALSE
                                     )
                                 ))
            if (!nrow(overlappingAnnotStart)) {
                stop(
                    'No overlap between CDS and transcripts were found. This is most likely due to a annoation problem around chromosome name.'
                )
            }

            # Annoate overlap ids
            overlappingAnnotStart$transcript_id <-
                localExons$transcript_id[overlappingAnnotStart$queryHits]
            overlappingAnnotStart$exon_id <- localExons$exon_id[
                overlappingAnnotStart$queryHits
                ]

            overlappingAnnotStart$cdsTranscriptID <- annotatedORFGR$id[
                overlappingAnnotStart$subjectHits
                ]
            overlappingAnnotStart$orf_id <- annotatedORFGR$orf_id[
                overlappingAnnotStart$subjectHits
                ]

            # subset to annoateted overlap
            overlappingAnnotStart <-
                overlappingAnnotStart[which(
                    overlappingAnnotStart$transcript_id ==
                        overlappingAnnotStart$cdsTranscriptID
                ), c('transcript_id',
                     'exon_id',
                     'cdsTranscriptID',
                     'orf_id')]

            # annoate with genomic site
            overlappingAnnotStart$orfGenomic <-
                start(annotatedORFGR)[match(
                    overlappingAnnotStart$orf_id, annotatedORFGR$orf_id
                )]


            ### Enrich exon information
            myExons <-
                as.data.frame(localExons[which(
                    localExons$transcript_id %in%
                        overlappingAnnotStart$transcript_id),])

            # Strand
            myExonPlus <- myExons[which(myExons$strand == '+'), ]
            myExonMinus <- myExons[which(myExons$strand == '-'), ]

            plusSplit <-
                split(myExonPlus$width, myExonPlus$transcript_id)
            minusSplit <-
                split(myExonMinus$width, myExonMinus$transcript_id)

            # cumsum
            myExonPlus$cumSum <-
                unlist(sapply(plusSplit , function(aVec) {
                    cumsum(c(0, aVec))[1:(length(aVec))]
                }))
            myExonMinus$cumSum <-
                unlist(sapply(minusSplit, function(aVec) {
                    cumsum(c(0, rev(aVec)))[(length(aVec)):1] # reverse
                }))

            # exon number
            myExonPlus$nrExon <-
                unlist(sapply(plusSplit, function(aVec) {
                    1:length(aVec)
                }))
            myExonMinus$nrExon <-
                unlist(sapply(minusSplit, function(aVec) {
                    1:length(aVec)
                }))

            # total nr exons
            myExonPlus$lastExonIndex <-
                unlist(sapply(plusSplit, function(aVec) {
                    rep(length(aVec), length(aVec))
                }))
            myExonMinus$lastExonIndex <-
                unlist(sapply(minusSplit, function(aVec) {
                    rep(1, length(aVec))
                }))

            # final exon exon junction trancipt position
            myExonPlus$finalJunctionPos <-
                unlist(sapply(plusSplit , function(aVec) {
                    rep(cumsum(c(0, aVec))[length(aVec)], times = length(aVec))
                }))
            myExonMinus$finalJunctionPos <-
                unlist(sapply(minusSplit, function(aVec) {
                    rep(cumsum(c(0, rev(
                        aVec
                    )))[length(aVec)], times = length(aVec))
                }))

            myExons2 <- rbind(myExonPlus, myExonMinus)

            ### Annoate with exon information
            matchIndex <-
                match(overlappingAnnotStart$exon_id, myExons2$exon_id)
            overlappingAnnotStart$strand <- myExons2$strand[matchIndex]
            overlappingAnnotStart$exon_start <- myExons2$start[matchIndex]
            overlappingAnnotStart$exon_end <- myExons2$end[matchIndex]
            overlappingAnnotStart$exon_cumsum <- myExons2$cumSum[matchIndex]
            overlappingAnnotStart$exon_nr <- myExons2$nrExon[matchIndex]
            overlappingAnnotStart$lastExonIndex <-
                myExons2$lastExonIndex[matchIndex]
            overlappingAnnotStart$finalJunctionPos <-
                myExons2$finalJunctionPos[matchIndex]

            ### Annoate with transcript coordinats
            overlappingAnnotStartPlus <-
                overlappingAnnotStart[which(
                    overlappingAnnotStart$strand == '+'), ]
            overlappingAnnotStartPlus$orfTranscript <-
                overlappingAnnotStartPlus$exon_cumsum + (
                    overlappingAnnotStartPlus$orfGenomic -
                        overlappingAnnotStartPlus$exon_start
                ) + 1
            overlappingAnnotStartPlus$junctionDistance <-
                overlappingAnnotStartPlus$finalJunctionPos -
                    overlappingAnnotStartPlus$orfTranscript + 3 # +3 because the ORF does not include the stop codon - but it should in this calculation

            overlappingAnnotStartMinus <-
                overlappingAnnotStart[which(
                    overlappingAnnotStart$strand == '-'), ]
            overlappingAnnotStartMinus$orfTranscript <-
                overlappingAnnotStartMinus$exon_cumsum + (
                    overlappingAnnotStartMinus$exon_end -
                        overlappingAnnotStartMinus$orfGenomic
                ) + 1
            overlappingAnnotStartMinus$junctionDistance <-
                overlappingAnnotStartMinus$finalJunctionPos -
                overlappingAnnotStartMinus$orfTranscript + 3 # +3 because the ORF does not include the stop codon - but it should in this calculation

            overlappingAnnotStart2 <-
                rbind(overlappingAnnotStartPlus,
                      overlappingAnnotStartMinus)
            overlappingAnnotStart2 <-
                overlappingAnnotStart2[order(
                    overlappingAnnotStart2$transcript_id,
                    overlappingAnnotStart2$exon_start,
                    overlappingAnnotStart2$exon_end
                ), ]

            ### devide into start and stop
            starInfo <-
                overlappingAnnotStart2[which(
                    grepl('^cds', overlappingAnnotStart2$orf_id)), ]
            stopInfo <-
                overlappingAnnotStart2[which(
                    grepl('^stop', overlappingAnnotStart2$orf_id)), ]

            ### predict PTC
            stopInfo$PTC <-
                stopInfo$exon_nr != stopInfo$lastExonIndex &
                stopInfo$junctionDistance > PTCDistance

            ### Merge the data
            starInfo2 <-
                unique(starInfo[, c('transcript_id',
                                    'orfGenomic',
                                    'exon_nr',
                                    'orfTranscript')])
            colnames(starInfo2) <-
                c('isoform_id',
                  'orfStartGenomic',
                  'orfStarExon',
                  'orfTransciptStart')

            stopInfo2 <-
                unique(stopInfo[, c(
                    'transcript_id',
                    'orfGenomic',
                    'exon_nr',
                    'orfTranscript',
                    'junctionDistance',
                    'lastExonIndex',
                    'PTC'
                )])
            colnames(stopInfo2) <-
                c(
                    'isoform_id',
                    'orfEndGenomic',
                    'orfEndExon',
                    'orfTransciptEnd',
                    'stopDistanceToLastJunction',
                    'stopIndex',
                    'PTC'
                )

            orfInfo <- merge(starInfo2, stopInfo2, by = 'isoform_id')
            orfInfo$orfTransciptLength  <-
                orfInfo$orfTransciptEnd - orfInfo$orfTransciptStart + 1

            # reorder
            orfInfo <-
                orfInfo[, c(
                    'isoform_id',
                    'orfTransciptStart',
                    'orfTransciptEnd',
                    'orfTransciptLength',
                    'orfStarExon',
                    'orfEndExon',
                    'orfStartGenomic',
                    'orfEndGenomic',
                    'stopDistanceToLastJunction',
                    'stopIndex',
                    'PTC'
                )]

            # make sure all ORFs are annotated (with NAs)
            orfInfo <-
                merge(orfInfo,
                      unique(myIsoAnot[, 'isoform_id', drop = FALSE]),
                      by = 'isoform_id',
                      all = TRUE)

        } else {
            # if no CDS were found
            warning(paste(
                'No CDS was found in the GTF file. Please make sure the GTF',
                'file have the CDS "feature" annotation. Adding NAs instead'
            ))

            orfInfo <- data.frame(
                isoform_id = unique(myIsoAnot$isoform_id),
                orfTransciptStart = NA,
                orfTransciptEnd = NA,
                orfTransciptLength = NA,
                orfStarExon = NA,
                orfEndExon = NA,
                orfStartGenomic = NA,
                orfEndGenomic = NA,
                stopDistanceToLastJunction = NA,
                stopIndex = NA,
                PTC = NA,
                stringsAsFactors = FALSE
            )
        }

        ### add to iso annotation
        myIsoAnot$PTC <-
            orfInfo$PTC[match(myIsoAnot$isoform_id, orfInfo$isoform_id)]

    }

    # Create exon_features grange
    myExons <-
        sort(mfGTF[exonAnoationIndex , c('transcript_id', 'gene_id')])
    colnames(myExons@elementMetadata) <- c('isoform_id', 'gene_id')

    # create replicates
    nrRep <-
        data.frame(
            condition = c('plaseholder1', 'plaseholder2'),
            nrReplicates = c(NA, NA),
            row.names = NULL,
            stringsAsFactors = FALSE
        )

    # create dummy feature
    repExp <- data.frame(
        isoform_id = myIsoAnot$isoform_id,
        plaseholder1 = NA,
        plaseholder2 = NA
    )

    designMatrix <-
        data.frame(
            sampleID = c('plaseholder1', 'plaseholder2'),
            condition = c('plaseholder1', 'plaseholder2')
        )

    ### Create switchList
    localSwichList <- createSwitchAnalyzeRlist(
        isoformFeatures = myIsoAnot,
        exons = myExons,
        designMatrix = designMatrix,
        isoformCountMatrix = repExp,
        sourceId = 'gtf'
    )

    if (addAnnotatedORFs) {
        # subset to those in list
        orfInfo <-
            orfInfo[which(orfInfo$isoform_id %in%
                              localSwichList$isoformFeatures$isoform_id), ]

        # check for negative ORF lengths
        isoformsToRemove <-
            orfInfo$isoform_id[which(orfInfo$orfTransciptLength < 0)]
        if (length(isoformsToRemove)) {
            genesToRemove <-
                localSwichList$isoformFeatures$gene_id[which(
                    localSwichList$isoformFeatures$isoform_id %in%
                        isoformsToRemove)]
            localSwichList <-
                subsetSwitchAnalyzeRlist(
                    localSwichList,
                    !localSwichList$isoformFeatures$gene_id %in% genesToRemove
                )

            warning(
                paste(
                    length(genesToRemove),
                    'genes where removed due to negative ORF lengths. This',
                    'typically occures because gene_id are not unique',
                    '(meaning are found multiple places accorss the genome).',
                    'Please note there might still be duplicated gene_id',
                    'located on the same chromosome.',
                    sep = ' '
                )
            )
        }
        localSwichList$orfAnalysis <- orfInfo[which(
            orfInfo$isoform_id %in% localSwichList$isoformFeatures$isoform_id)
            ,]
    }

    # Return SpliceRList
    return(localSwichList)
}

importIsoformExpression <- function(
    parentDir,
    calculateCountsFromAbundance=TRUE,
    interLibNormTxPM=TRUE,
    normalizationMethod='TMM',
    pattern='',
    invertPattern=FALSE,
    ignore.case=FALSE,
    showProgress = TRUE,
    quiet = FALSE
) {
    ### To do
    # Could get a "summarize to gene level" option

    ### Initialize
    if(TRUE) {
        if( !normalizationMethod %in% c("TMM", "RLE", "upperquartile") ){
            stop('Metod supplied to "normalizationMethod" must be one of "TMM", "RLE" or "upperquartile". See documentation of edgeR::calcNormFactors for more info')
        }

        analysisCount <- 2 +
            as.integer(interLibNormTxPM)

        if (showProgress &  !quiet) {
            progressBar <- 'text'
            progressBarLogic <- TRUE
        } else {
            progressBar <- 'none'
            progressBarLogic <- FALSE
        }

        ### data.frame with nesseary info
        supportedTypes <- data.frame(
            orign          = c('Kallisto'       , 'Salmon'         , 'RSEM'),
            fileName       = c('abundance.tsv'  , 'quant.sf'       , 'isoforms.results'),
            id             = c('target_id'      , 'Name'           , 'transcript_id'),
            countCol       = c('est_counts'     , 'NumReads'       , 'expected_count'),
            tpmCol         = c('tpm'            , 'TPM'            , 'TPM'),
            eLengthCol     = c('eff_length'     , 'EffectiveLength', 'effective_length'),
            stringsAsFactors = FALSE
        )
    }

    ### Identify directories of interest
    if (TRUE) {
        if (!quiet) {
            message('Step 1 of ', analysisCount, ': Identifying which algorithm was used...')
        }
        dirList <- split(
            list.dirs(
                path = parentDir,
                full.names = FALSE,
                recursive = FALSE
            ),
            list.dirs(
                path = parentDir,
                full.names = FALSE,
                recursive = FALSE
            )
        )
        dirList <- dirList[which(sapply(dirList, nchar) > 0)]

        ### Extract those where there are files of interest
        dirList <-
            dirList[sapply(
                dirList,
                FUN = function(aDir) {
                    # aDir <- dirList[[1]]
                    localFiles <-
                        list.files(
                            paste0(parentDir, '/', aDir),
                            recursive = FALSE
                        )

                    any(sapply(
                        paste(supportedTypes$fileName, '$', sep = ''),
                        function(aFileName) {
                            grepl(pattern = aFileName, x =  localFiles)
                        }))
                }
            )]

        if (length(dirList) == 0) {
            stop('No directories of interest were found')
        }

    }

    ### Identify input type
    if(TRUE) {
        dataAnalyed <- supportedTypes[which(
            sapply(
                paste0(supportedTypes$fileName,'$'),
                function(aFileName) {
                    any(grepl(
                        pattern = aFileName,
                        x = list.files(paste0( parentDir, '/', dirList[[1]] ))
                    ))
                })
        ), ]

        if (nrow(dataAnalyed) > 1) {
            stop('Could not uniquely identify file type - please contact developer')
        }
        if (!quiet) {
            message(paste('    The quantification algorithm used was:', dataAnalyed$orign, sep = ' '))
        }

    }

    ### Import files with txtimport
    if(TRUE) {
        if (!quiet) {
            message('Step 2 of ', analysisCount, ': Reading data...')
        }

        # make vector with paths
        localFiles <- paste0(
            parentDir,
            '/',
            unlist(dirList),
            '/',
            dataAnalyed$fileName
        )
        names(localFiles) <- names(dirList)

        ### Subset to those of interest
        if( invertPattern ) {
            localFiles <- localFiles[which(
                ! grepl(
                    pattern = pattern,
                    x = localFiles,
                    ignore.case=ignore.case
                )
            )]
        } else {
            localFiles <- localFiles[which(
                grepl(
                    pattern = pattern,
                    x = localFiles,
                    ignore.case=ignore.case
                )
            )]
        }

        ### Use Txtimporter to import data
        if (!quiet) {
            localDataList <- tximport(
                files = localFiles,
                type = tolower(dataAnalyed$orign),
                txOut = TRUE, # to get isoform expression
                countsFromAbundance = ifelse(
                    test = calculateCountsFromAbundance,
                    yes= 'lengthScaledTPM',
                    no='no'
                )
            )
        } else {
            suppressMessages(
                localDataList <- tximport(
                    files = localFiles,
                    type = tolower(dataAnalyed$orign),
                    txOut = TRUE, # to get isoform expression
                    countsFromAbundance = ifelse(
                        test = calculateCountsFromAbundance,
                        yes= 'lengthScaledTPM',
                        no='no'
                    )
                )
            )
        }

    }

    analysisDone <- 2

    ### Noralize TxPM values based on effective counts
    if(interLibNormTxPM) {
        if (!quiet) {
            message('Step ', analysisDone, ' of ', analysisCount, ': Normalizing TxPM values via edgeR...')
        }

        ### calclate new coints
        newCounts <- localDataList$abundance * localDataList$length

        ### Scale new to same total counts as org matrix
        countsSum <- colSums(localDataList$counts)
        newSum <- colSums(newCounts)
        countsMat <- t(t(newCounts) * (countsSum/newSum))

        ### Calculate normalization factors
        localDGE <- edgeR::DGEList(countsMat, remove.zeros = TRUE)
        localDGE <- edgeR::calcNormFactors(localDGE, method = normalizationMethod)

        ### Apply normalization factors
        localDataList$abundance <- t(t(localDataList$abundance) / localDGE$samples$norm.factors)
    }

    ### Final masssageing
    if(TRUE) {
        ### massage
        localDataList$abundance <- as.data.frame(localDataList$abundance)
        localDataList$counts <- as.data.frame(localDataList$counts)
        localDataList$length <- as.data.frame(localDataList$length)

        localDataList$countsFromAbundance <- NULL

        reorderCols <- function(x) {
            x[,c( ncol(x), 1:(ncol(x)-1) )]
        }

        localDataList$abundance <- reorderCols( localDataList$abundance)
        localDataList$counts    <- reorderCols( localDataList$counts   )
        localDataList$length    <- reorderCols( localDataList$length   )

        ### Add options
        localDataList$importOptions <- list(
            'calculateCountsFromAbundance'= calculateCountsFromAbundance,
            'interLibNormTxPM'= interLibNormTxPM,
            'normalizationMethod'= normalizationMethod
        )

        if (!quiet) {
            message('Done\n')
        }
    }

    return(localDataList)
}

importRdata <- function(
    isoformCountMatrix,
    isoformRepExpression = NULL,
    designMatrix,
    isoformExonAnnoation,
    comparisonsToMake = NULL,
    addAnnotatedORFs = FALSE,
    onlyConsiderFullORF = FALSE,
    removeNonConvensionalChr = FALSE,
    PTCDistance = 50,
    foldChangePseudoCount = 0.01,
    showProgress = TRUE,
    quiet = FALSE
) {
    ### Set up progress
    if (showProgress &  !quiet) {
        progressBar <- 'text'
        progressBarLogic <- TRUE
    } else {
        progressBar <- 'none'
        progressBarLogic <- FALSE
    }


    ### Test whether imput data fits together
    if (TRUE) {
        ### Contains the collums they should
        if (TRUE) {
            ### Colnames
            if (!any(colnames(isoformCountMatrix) == 'isoform_id')) {
                stop(paste(
                    'The data.frame passed to the \'isoformCountMatrix\'',
                    'argument must contain a \'isoform_id\' column'
                ))
            }

            if (!is.null(isoformRepExpression)) {
                if (!any(colnames(isoformRepExpression) == 'isoform_id')) {
                    stop(paste(
                        'The data.frame passed to the \'isoformCountMatrix\'',
                        'argument must contain a \'isoform_id\' column'
                    ))
                }
            }

            if (!all(c('sampleID', 'condition') %in% colnames(designMatrix))) {
                stop(paste(
                    'The data.frame passed to the \'isoformRepExpression\'',
                    'argument must contain both a \'sampleID\' and a',
                    '\'condition\' column'
                ))
            }
            if (length(unique(designMatrix$condition)) < 2) {
                stop('The supplied \'designMatrix\' only contains 1 condition')
            }
            if (!is.null(comparisonsToMake)) {
                if (!all(c('condition_1', 'condition_2') %in%
                         colnames(comparisonsToMake))) {
                    stop(paste(
                        'The data.frame passed to the \'comparisonsToMake\'',
                        'argument must contain both a \'condition_1\' and a',
                        '\'condition_2\' column indicating',
                        'the comparisons to make'
                    ))
                }
            }
        }

        ### Removel potential levels
        if (TRUE) {
            designMatrix$condition <- as.character(designMatrix$condition)
            designMatrix$sampleID  <-
                as.character(designMatrix$sampleID)

            if (!is.null(comparisonsToMake)) {
                comparisonsToMake$condition_1 <-
                    as.character(comparisonsToMake$condition_1)
                comparisonsToMake$condition_2 <-
                    as.character(comparisonsToMake$condition_2)
            }

            if (!is.null(isoformRepExpression)) {
                isoformRepExpression$isoform_id <-
                    as.character(isoformRepExpression$isoform_id)
            }
            isoformCountMatrix$isoform_id <-
                as.character(isoformCountMatrix$isoform_id)
        }

        ### Check if it fits togehter
        if (TRUE) {
            if (!all(designMatrix$sampleID %in% colnames(isoformCountMatrix))) {
                stop(paste(
                    'Each sample stored in \'designMatrix$sampleID\' must have',
                    'a corresponding expression column in \'isoformCountMatrix\''
                ))
            }
            if (!is.null(isoformRepExpression)) {
                if (!all(designMatrix$sampleID %in%
                         colnames(isoformRepExpression))) {
                    stop(paste(
                        'Each sample stored in \'designMatrix$sampleID\' must',
                        'have a corresponding expression column',
                        'in \'isoformRepExpression\''
                    ))
                }
            }
            if (!all(designMatrix$sampleID %in% colnames(isoformCountMatrix))) {
                stop(
                    'Each sample stored in \'designMatrix$sampleID\' must have',
                    'a corresponding expression column in \'isoformRepExpression\''
                )
            }

            if (!is.null(comparisonsToMake)) {
                if (!all(
                    c(
                        comparisonsToMake$condition_1,
                        comparisonsToMake$condition_2
                    ) %in% designMatrix$condition
                )) {
                    stop(paste(
                        'The conditions supplied in comparisonsToMake and',
                        'designMatrix does not match'
                    ))
                }
            } else {
                # create it
                comparisonsToMake <-
                    allPairwiseFeatures(designMatrix$condition)
                colnames(comparisonsToMake) <-
                    c('condition_1', 'condition_2')
            }
        }
    }

    ### Handle annotation imput
    if (!quiet) {
        message('Step 1 of 5: Obtaining annotation...')
    }
    if (TRUE) {
        ### Import GTF is nessesary
        if (class(isoformExonAnnoation) == 'character') {
            gtfImported <- TRUE

            if (length(isoformExonAnnoation) != 1) {
                stop(paste(
                    'The path supplied to isoformExonAnnoation',
                    'must be of length 1'
                ))
            }
            if (!grepl('gtf$', isoformExonAnnoation)) {
                warning(paste(
                    'The GTF file supplied to isoformExonAnnoation',
                    'does not contain the suffix \'.gtf\''
                ))
            }
            if (!file.exists(isoformExonAnnoation)) {
                stop(paste(
                    'The file paht supplied to isoformExonAnnoation',
                    'points to a file that does not exists'
                ))
            }

            if (!quiet) {
                message('importing GTF (this may take a while)')
            }
            gtfSwichList <- importGTF(
                pathToGTF = isoformExonAnnoation,
                addAnnotatedORFs = addAnnotatedORFs,
                onlyConsiderFullORF = onlyConsiderFullORF,
                PTCDistance = PTCDistance,
                removeNonConvensionalChr = removeNonConvensionalChr,
                quiet = TRUE
            )

            ### Extract wanted annotation files form the spliceR object
            isoformExonStructure <-
                gtfSwichList$exons[, c('isoform_id', 'gene_id')]
            isoformExonStructure <- sort(isoformExonStructure)

            isoformAnnotation <-
                unique(gtfSwichList$isoformFeatures[,c(
                    'isoform_id', 'gene_id', 'gene_name'
                )])

            if (addAnnotatedORFs & gtfImported) {
                isoORF <- gtfSwichList$orfAnalysis
            }
        } else {
            gtfImported <- FALSE

            ### Devide the data
            isoformExonStructure <-
                isoformExonAnnoation[, c('isoform_id', 'gene_id')]

            isoformAnnotation <-
                unique(as.data.frame(isoformExonAnnoation@elementMetadata))
            if (!'gene_name' %in% colnames(isoformAnnotation)) {
                isoformAnnotation$gene_name <- NA
            }
        }

        ### Reduce to those expressed
        isoformExonStructure <- isoformExonStructure[which(
            isoformExonStructure$isoform_id %in% isoformCountMatrix$isoform_id
        ), ]
        isoformAnnotation <-isoformAnnotation[which(
            isoformAnnotation$isoform_id    %in% isoformCountMatrix$isoform_id
        ), ]

        if (length(isoformExonStructure) == 0 |
            nrow(isoformAnnotation) == 0) {
            stop('The annotation does not fit the expression data')
        }

        ### Test the obtained annoation
        if (!all(c('isoform_id', 'gene_id', 'gene_name') %in%
                 colnames(isoformAnnotation))) {
            stop(paste(
                'The data.frame passed to the \'isoformAnnotation\' argument',
                'must contain the following columns \'isoform_id\',',
                '\'gene_id\' and \'gene_name\''
            ))
        }
        if (any(is.na(isoformAnnotation[, c('isoform_id', 'gene_id')]))) {
            stop(paste(
                'The \'isoform_id\' and \'gene_id\' columns in the data.frame',
                'passed to the \'isoformAnnotation\' argument are not allowed',
                'to contain NAs'
            ))
        }
        if (!'isoform_id' %in% colnames(isoformExonStructure@elementMetadata)) {
            stop(paste(
                'The GenomicRanges (GRanges) object passed to the',
                '\'isoformExonStructure\' argument must contain both a',
                '\'isoform_id\' and \'gene_id\' metadata column'
            ))
        }

        if (jaccardDistance(isoformCountMatrix$isoform_id,
                            isoformAnnotation$isoform_id) == 0) {
            stop('The annotation does not fit the expression data')
        }
        if (jaccardDistance(isoformCountMatrix$isoform_id,
                            isoformExonStructure$isoform_id) == 0) {
            stop('The annotation does not fit the expression data')
        }
    }

    ### Reduce to used data
    if (TRUE) {
        designMatrix <-
            designMatrix[which(
                designMatrix$condition %in% c(
                    comparisonsToMake$condition_1,
                    comparisonsToMake$condition_2
                )
            ), ]
        designMatrix$sampleID <- as.character(designMatrix$sampleID)
        designMatrix$condition <-
            as.character(designMatrix$condition)

        if (!is.null(isoformRepExpression)) {
            isoformRepExpression <-
                isoformRepExpression[, which(
                    colnames(isoformRepExpression) %in%
                        c('isoform_id', designMatrix$sampleID)
                )]
        }
        isoformCountMatrix <-
            isoformCountMatrix[, which(
                colnames(isoformCountMatrix) %in%
                    c('isoform_id', designMatrix$sampleID))]
        isoformCountMatrix <-
            isoformCountMatrix[,c(
                which(colnames(isoformCountMatrix) == 'isoform_id'),
                which(colnames(isoformCountMatrix) != 'isoform_id')
            )]
        rownames(isoformCountMatrix) <- NULL
    }

    ### If nessesary calculate RPKM values
    if (!quiet) {
        message('Step 2 of 5: Calculating gene expression...')
    }
    if (TRUE) {
        if (is.null(isoformRepExpression)) {
            ### Extract isoform lengths
            isoformLengths <- sapply(
                X = split(
                    isoformExonStructure@ranges@width,
                    f = isoformExonStructure$isoform_id
                ),
                FUN = sum
            )

            ### Calulate TPM
            # convert to matrix
            localCM <- isoformCountMatrix
            rownames(localCM) <- localCM$isoform_id
            localCM$isoform_id <- NULL
            localCM <- as.matrix(localCM)

            myTPM <- t(t(localCM) / colSums(localCM)) * 1e6

            ### Calculate RPKM
            isoformLengths <-
                isoformLengths[match(rownames(myTPM), names(isoformLengths))]

            isoformRepExpression <-
                as.data.frame(myTPM / (isoformLengths / 1e3))
            isoformRepExpression$isoform_id <-
                rownames(isoformRepExpression)
            isoformRepExpression <-
                isoformRepExpression[, c(
                    which(colnames(isoformRepExpression) == 'isoform_id'),
                    which(colnames(isoformRepExpression) != 'isoform_id')
                )]
            rownames(isoformRepExpression) <- NULL
        }
    }

    ### Remove isoforms not annoated
    if (TRUE) {
        diffAnnot <-
            setdiff(isoformRepExpression$isoform_id ,
                    isoformAnnotation$isoform_id)
        if (length(diffAnnot)) {
            warning(
                paste(
                    'There were',
                    length(diffAnnot),
                    'isoforms whoes expression were measured',
                    'but which were not annotated.',
                    'These are included in the gene expression estimation',
                    'but not in the final switchAnalyzeRlist',
                    sep = ' '
                )
            )
            isoformRepExpression <-
                isoformRepExpression[which(
                    isoformRepExpression$isoform_id %in%
                        isoformAnnotation$isoform_id), ]
        }
    }

    ### Sum to gene level gene expression - updated
    if(TRUE) {
        ### add gene_id
        isoformRepExpression2 <- isoformRepExpression
        isoformRepExpression2$gene_id <-
            isoformAnnotation$gene_id[match(isoformRepExpression2$isoform_id,
                                            isoformAnnotation$isoform_id)]

        ### Sum to gene level
        geneRepExpression <- isoformToGeneExp(
            isoformRepExpression2,
            showProgress = showProgress,
            quiet = quiet
        )
    }

    ### in each condition analyzed get mean and standard error of gene and isoforms
    if (!quiet) {
        message('Step 3 of 5: Merging gene and isoform expression...')
    }
    if (TRUE) {
        conditionList <-
            split(designMatrix$sampleID, f = designMatrix$condition)
        conditionSummary <-
            llply(
                .data = conditionList,
                .progress = progressBar,
                .fun = function(sampleVec) {
                    # sampleVec <- conditionList[[1]]
                    ### Isoform
                    isoIndex <-
                        which(colnames(isoformRepExpression) %in% sampleVec)

                    isoSummary <- data.frame(
                        isoform_id = isoformRepExpression$isoform_id,
                        iso_value = rowMeans(isoformRepExpression[, isoIndex]),
                        iso_std = apply(isoformRepExpression[, isoIndex], 1, sd),
                        stringsAsFactors = FALSE
                    )
                    isoSummary$iso_stderr <-
                        isoSummary$iso_std / sqrt(length(sampleVec))
                    isoSummary$iso_std <- NULL

                    ### Gene
                    geneIndex <-
                        which(colnames(geneRepExpression) %in% sampleVec)

                    geneSummary <- data.frame(
                        gene_id = geneRepExpression$gene_id,
                        gene_value = rowMeans(geneRepExpression[, geneIndex]),
                        gene_std = apply(geneRepExpression[, geneIndex], 1, sd),
                        stringsAsFactors = FALSE
                    )
                    geneSummary$gene_stderr <-
                        geneSummary$gene_std / sqrt(length(sampleVec))
                    geneSummary$gene_std <- NULL

                    ### Combine
                    combinedData <-
                        merge(isoformAnnotation, geneSummary, by = 'gene_id')
                    combinedData <-
                        merge(combinedData, isoSummary, by = 'isoform_id')
                    ### return result
                    return(combinedData)
                }
            )
    }

    ### Use comparisonsToMake to create the isoform comparisons
    if (!quiet) {
        message('Step 4 of 5: Making comparisons...')
    }
    if (TRUE) {
        isoAnnot <-
            ddply(
                .data = comparisonsToMake,
                .variables = c('condition_1', 'condition_2'),
                .drop = TRUE,
                .progress = progressBar,
                .fun = function(aDF) {
                    ### Extract data
                    cond1data <- conditionSummary[[aDF$condition_1]]
                    cond2data <- conditionSummary[[aDF$condition_2]]

                    # modify colnames
                    matchIndex <-
                        match(
                            c(
                                'gene_value',
                                'gene_stderr',
                                'iso_value',
                                'iso_stderr'
                            ),
                            colnames(cond1data)
                        )
                    colnames(cond1data)[matchIndex] <-
                        paste(colnames(cond1data)[matchIndex], '_1', sep = '')

                    matchIndex <-
                        match(
                            c(
                                'gene_value',
                                'gene_stderr',
                                'iso_value',
                                'iso_stderr'
                            ),
                            colnames(cond2data)
                        )
                    colnames(cond2data)[matchIndex] <-
                        paste(colnames(cond2data)[matchIndex], '_2', sep = '')

                    combinedIso <- merge(cond1data,
                                         cond2data[, c(
                                             'isoform_id',
                                             'gene_value_2',
                                             'gene_stderr_2',
                                             'iso_value_2',
                                             'iso_stderr_2'
                                         )],
                                         by = 'isoform_id')

                    ### Add comparison data
                    combinedIso$condition_1 <- aDF$condition_1
                    combinedIso$condition_2 <- aDF$condition_2
                    return(combinedIso)
                }
            )

        ### Add comparison data
        # Log2FC
        ps <- foldChangePseudoCount

        isoAnnot$gene_log2_fold_change <-
            log2((isoAnnot$gene_value_2 + ps) / (isoAnnot$gene_value_1 + ps))
        isoAnnot$iso_log2_fold_change  <-
            log2((isoAnnot$iso_value_2  + ps) / (isoAnnot$iso_value_1  + ps))

        # qValues
        isoAnnot$gene_q_value <- NA
        isoAnnot$iso_q_value  <- NA

        # Isoform fraction values
        isoAnnot$IF1 <-
            round(isoAnnot$iso_value_1 / isoAnnot$gene_value_1 , digits = 4)
        isoAnnot$IF2 <-
            round(isoAnnot$iso_value_2 / isoAnnot$gene_value_2 , digits = 4)
        isoAnnot$dIF <- isoAnnot$IF2 - isoAnnot$IF1

        # Swich values
        isoAnnot$isoform_switch_q_value <- NA
        isoAnnot$gene_switch_q_value    <- NA

        ### Sort
        matchVector <-
            c(
                'isoform_id',
                'gene_id',
                'condition_1',
                'condition_2',
                'gene_name',
                'class_code',
                'gene_value_1',
                'gene_value_2',
                'gene_stderr_1',
                'gene_stderr_2',
                'gene_log2_fold_change',
                'gene_q_value',
                'iso_value_1',
                'iso_value_2',
                'iso_stderr_1',
                'iso_stderr_2',
                'iso_log2_fold_change',
                'iso_q_value',
                'IF1',
                'IF2',
                'dIF',
                'isoform_switch_q_value',
                'gene_switch_q_value'
            )
        matchVector <-
            na.omit(match(matchVector, colnames(isoAnnot)))

        isoAnnot <- isoAnnot[, matchVector]
    }

    ### Create the swichList
    if (!quiet) {
        message('Step 5 of 5: Making switchAnalyzeRlist object...')
    }
    if (TRUE) {
        ### Create switchList
        dfSwichList <- createSwitchAnalyzeRlist(
            isoformFeatures = isoAnnot,
            exons = isoformExonStructure,
            designMatrix = designMatrix,
            isoformCountMatrix = isoformCountMatrix,
            sourceId = 'data.frames'
        )

        ### Add orf if extracted
        if (addAnnotatedORFs & gtfImported) {
            dfSwichList$isoformFeatures$PTC <-
                isoORF$PTC[match(dfSwichList$isoformFeatures$isoform_id,
                                 isoORF$isoform_id)]

            isoORF <-
                isoORF[which(isoORF$isoform_id %in%
                                 isoformRepExpression$isoform_id), ]

            dfSwichList$orfAnalysis <- isoORF
        }
    }

    if (!quiet) {
        message('Done')
    }
    return(dfSwichList)
}


### Prefilter
preFilter <- function(
    switchAnalyzeRlist,
    geneExpressionCutoff = 1,
    isoformExpressionCutoff = 0,
    IFcutoff = 0.01,
    geneCutoffInBothConditions = TRUE,
    acceptedIsoformClassCode = NULL,
    removeSingleIsoformGenes = TRUE,
    reduceToSwitchingGenes = FALSE,
    keepIsoformInAllConditions = FALSE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    quiet = FALSE
) {
    ### Test input
    if (TRUE) {
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(paste(
                'The object supplied to \'switchAnalyzeRlist\'',
                'must be a \'switchAnalyzeRlist\''
            ))
        }

        if (!is.null(isoformExpressionCutoff)) {
            if (!is.numeric(isoformExpressionCutoff)) {
                stop('The isoformExpressionCutoff argument must be a numeric')
            }
        }
        if (!is.null(geneExpressionCutoff)) {
            if (!is.numeric(geneExpressionCutoff)) {
                stop('The geneExpressionCutoff argument must be a numeric')
            }
        }
        if (!is.null(IFcutoff)) {
            if (!is.numeric(IFcutoff)) {
                stop('The IFcutoff argument must be a numeric')
            }
        }

        if (!is.null(acceptedIsoformClassCode)) {
            if (!'class_code' %in%
                colnames(switchAnalyzeRlist$isoformFeatures)) {
                stop(paste(
                    'The filter on class codes can only be used if',
                    'the switchAnalyzeRlist was generated from cufflinks data'
                ))
            }
        }

        if (!is.logical(removeSingleIsoformGenes)) {
            stop('The removeSingleIsoformGenes must be either TRUE or FALSE')
        }
        if (!is.logical(geneCutoffInBothConditions)) {
            stop('The geneCutoffInBothConditions must be either TRUE or FALSE')
        }

        if (reduceToSwitchingGenes) {
            if (alpha < 0 |
                alpha > 1) {
                stop('The alpha parameter must be between 0 and 1 ([0,1]).')
            }
            if (alpha > 0.05) {
                warning(paste(
                    'Most journals and scientists consider an alpha larger',
                    'than 0.05 untrustworthy. We therefore recommend using',
                    'alpha values smaller than or queal to 0.05'
                ))
            }
        }
    }

    ### Find which isoforms to keep
    if (TRUE) {
        ### Extract data
        columnsToExtraxt <-
            c(
                'iso_ref',
                'gene_ref',
                'isoform_id',
                'gene_id',
                'class_code',
                'gene_value_1',
                'gene_value_2',
                'iso_value_1',
                'iso_value_2',
                'IF1',
                'IF2',
                'dIF',
                'gene_switch_q_value'
            )
        columnsToExtraxt <-
            na.omit(match(
                columnsToExtraxt,
                colnames(switchAnalyzeRlist$isoformFeature)
            ))
        #localData <- unique( switchAnalyzeRlist$isoformFeature[, columnsToExtraxt ] ) # no need
        localData <-
            switchAnalyzeRlist$isoformFeature[, columnsToExtraxt]

        # Count features
        transcriptCount <- length(unique(localData$isoform_id))

        ### Do filtering
        if (reduceToSwitchingGenes &
            any(!is.na(localData$gene_switch_q_value))) {
            localData <- localData[which(localData$gene_switch_q_value < alpha &
                                             abs(localData$dIF) > dIFcutoff), ]
            if (!nrow(localData)) {
                stop('No genes were left after filtering for switching genes')
            }
        }

        if (!is.null(geneExpressionCutoff)) {
            if (geneCutoffInBothConditions) {
                localData <- localData[which(
                    localData$gene_value_1 > geneExpressionCutoff &
                        localData$gene_value_2 > geneExpressionCutoff
                ), ]
            } else {
                localData <- localData[which(
                    localData$gene_value_1 > geneExpressionCutoff |
                        localData$gene_value_2 > geneExpressionCutoff
                ), ]
            }
            if (!nrow(localData)) {
                stop('No genes were left after filtering for gene expression')
            }
        }

        if (!is.null(isoformExpressionCutoff)) {
            localData <- localData[which(
                localData$iso_value_1 > isoformExpressionCutoff |
                    localData$iso_value_2 > isoformExpressionCutoff
            ), ]
            if (!nrow(localData)) {
                stop('No genes were left after filtering for isoform expression')
            }
        }

        if (!is.null(isoformExpressionCutoff)) {
            localData <- localData[which(localData$IF1 > IFcutoff |
                                             localData$IF2 > IFcutoff), ]
            if (!nrow(localData)) {
                stop('No genes were left after filtering for isoform fraction (IF) values')
            }
        }

        if (!is.null(acceptedIsoformClassCode)) {
            localData <- localData[which(localData$class_code %in% acceptedIsoformClassCode), ]

            if (!nrow(localData)) {
                stop('No genes were left after filtering for isoform class codes')
            }
        }

        if (removeSingleIsoformGenes) {
            transcriptsPrGene <-
                split(localData$iso_ref, f = localData$gene_ref)
            transcriptsPrGene <- lapply(transcriptsPrGene, unique)

            genesToKeep <-
                names(transcriptsPrGene)[which(sapply(transcriptsPrGene, function(x)
                    length(x) > 1))]

            if (!length(genesToKeep)) {
                stop('No genes were left after filtering for mutlipe transcrip genes')
            }

            localData <-
                localData[which(localData$gene_ref %in% genesToKeep),]
        }

    }

    ### Do filtering
    if (keepIsoformInAllConditions) {
        switchAnalyzeRlist <- subsetSwitchAnalyzeRlist(
            switchAnalyzeRlist,
            switchAnalyzeRlist$isoformFeatures$isoform_id %in%
                localData$isoform_id
        )
    } else {
        switchAnalyzeRlist <- subsetSwitchAnalyzeRlist(
            switchAnalyzeRlist,
            switchAnalyzeRlist$isoformFeatures$iso_ref %in% localData$iso_ref
        )
    }

    ### Repport filtering
    transcriptsLeft <- length(unique(localData$isoform_id))
    transcriptsRemoved <- transcriptCount - transcriptsLeft
    percentRemoved <-
        round(transcriptsRemoved / transcriptCount, digits = 4) * 100

    if (!quiet) {
        message(
            paste(
                'The fitering removed ',
                transcriptsRemoved,
                ' ( ',
                percentRemoved,
                '% of ) transcripts. There is now ',
                transcriptsLeft,
                ' isoforms left',
                sep = ''
            )
        )
    }

    ### return result
    return(switchAnalyzeRlist)
}
