# Helper functions
evalSig <- function(pValue, alphas) {
    sapply(pValue, function(x) {
        if( is.na(x) ) {
            return('NA')
        } else if( x < min( alphas) ) {
            sigLevel <- '***'
        } else if ( x < max( alphas) ) {
            sigLevel <- '*'
        } else {
            sigLevel <- 'ns'
        }
        return(sigLevel)
    })
}

medianQuartile <- function(x){
    out <- quantile(x, probs = c(0.25,0.5,0.75))
    names(out) <- c("ymin","y","ymax")
    return(out)
}

extractSigData <- function(
    switchAnalyzeRlist,
    alpha=0.05,
    dIFcutoff = 0.1,
    log2FCcutoff = 1,
    featureToExtract = 'isoformUsage'
) {
    ### Test input
    # done by the parrent functions
    if (featureToExtract == 'all') {
        return(switchAnalyzeRlist$isoformFeatures$iso_ref)
    }

    ### Extract sig iso usage
    if (featureToExtract == 'isoformUsage') {
        idsToExtract <-
            switchAnalyzeRlist$isoformFeatures$iso_ref[which(
                switchAnalyzeRlist$isoformFeatures$gene_switch_q_value < alpha &
                    abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
            )]

        return(idsToExtract)
    }
    if (featureToExtract == 'isoformExp') {
        idsToExtract <- switchAnalyzeRlist$isoformFeatures$iso_ref[which(
            switchAnalyzeRlist$isoformFeatures$iso_q_value < alpha &
                abs(
                    switchAnalyzeRlist$isoformFeatures$iso_log2_fold_change
                ) > log2FCcutoff
        )]
        return(idsToExtract)
    }
    if (featureToExtract == 'geneExp') {
        idsToExtract <- switchAnalyzeRlist$isoformFeatures$iso_ref[which(
            switchAnalyzeRlist$isoformFeatures$gene_q_value < alpha &
                abs(
                    switchAnalyzeRlist$isoformFeatures$gene_log2_fold_change
                ) > log2FCcutoff
        )]
        return(idsToExtract)
    }
}

# acutal function
extractGenomeWideAnalysis <- function(
    switchAnalyzeRlist,
    featureToExtract = 'isoformUsage',
    annotationToAnalyze = 'all',
    alpha=0.05,
    dIFcutoff = 0.1,
    log2FCcutoff = 1,
    violinPlot=TRUE,
    alphas=c(0.05, 0.001),
    localTheme=theme_bw(),
    plot=TRUE,
    returnResult=TRUE
) {
    ### Test input
    if (TRUE) {
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(
                'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
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
        if (dIFcutoff < 0 | dIFcutoff > 1) {
            stop('The dIFcutoff must be in the interval [0,1].')
        }
        if (log2FCcutoff < 0) {
            stop(
                'The log2FCcutoff cannot be negative (as the cutoff is applied to absolute values)'
            )
        }
        if (length(alphas) != 2) {
            stop('A vector of length 2 must be given to the argument \'alphas\'')
        }
        if (any(alphas < 0) |
            any(alphas > 1)) {
            stop('The \'alphas\' parameter must be numeric between 0 and 1 ([0,1]).')
        }
        if (any(alphas > 0.05)) {
            warning(
                'Most journals and scientists consider an alpha larger than 0.05 untrustworthy. We therefore recommend using alpha values smaller than or queal to 0.05'
            )
        }
        if (!featureToExtract %in%
            c('isoformUsage', 'isoformExp', 'geneExp', 'all')
        ) {
            stop(
                'The \'featureToExtract\' argument must be  \'isoformUsage\', \'isoformExp\',  \'geneExp\' or \'all\''
            )
        }

        ### Check annotation
        okAnnot <-
            c(
                'ORF',
                'NMD_status',
                'coding_potential',
                'signal_peptide_identified',
                'domains_identified',
                'intron_retention',
                'switch_consequences',
                'isoform_class_code'
            )
        if ('all' %in% annotationToAnalyze) {
            annotationToAnalyze <- okAnnot
        }
        if (!all(annotationToAnalyze %in% okAnnot)) {
            stop(
                paste(
                    'The \'annotationToAnalyze\' argument must be a one (or multiple) of: \'',
                    paste(annotationToAnalyze, collapse = '\', \''),
                    '\'',
                    sep = ''
                )
            )
        }

        ### Replace with annotaion names with collum names
        annotToExtract <-
            unique(gsub('NMD_status|ORF', 'PTC', annotationToAnalyze))
        annotToExtract <- gsub(
            'coding_potential',
            'codingPotential',
            annotToExtract
        )
        annotToExtract <- gsub(
            'domains_identified'  ,
            'domain_identified',
            annotToExtract
        )
        annotToExtract <-gsub(
            'intron_retention',
            'IR',
            annotToExtract
        )
        annotToExtract <- gsub(
            'switch_consequences' ,
            'switchConsequencesGene',
            annotToExtract
        )
        annotToExtract <- gsub(
            'isoform_class_code',
            'class_code',
            annotToExtract
        )

        if (length(annotToExtract) == 0) {
            stop(
                'Somthing in the annoation decoding went wrong - please send a small example dataset reconstructing the mistake to the developers.'
            )
        }
    }

    ### Extract annotation
    if (TRUE) {
        switchAnalyzeRlist$isoformFeatures$comparison <- paste(
            switchAnalyzeRlist$isoformFeatures$condition_1,
            switchAnalyzeRlist$isoformFeatures$condition_2,
            sep = ' vs '
        )

        isoformsToAnalyze <- extractSigData(
            switchAnalyzeRlist = switchAnalyzeRlist,
            alpha = alpha,
            dIFcutoff = dIFcutoff,
            log2FCcutoff = log2FCcutoff,
            featureToExtract = featureToExtract
        )

        colToExtract <-
            c('iso_ref',
              'isoform_id',
              'comparison',
              'IF1',
              'IF2',
              annotToExtract)
        colToExtract <-
            intersect(colToExtract,
                      colnames(switchAnalyzeRlist$isoformFeatures))

        isoData <- switchAnalyzeRlist$isoformFeatures[
            which(
                switchAnalyzeRlist$isoformFeatures$iso_ref %in%
                    isoformsToAnalyze
            ),
            na.omit(match(
                colToExtract ,
                colnames(switchAnalyzeRlist$isoformFeatures)
            ))
        ]



        if (nrow(isoData) == 0) {
            stop('No data left after filtering')
        }
    }

    ### Overwrite annotation with propper categories
    if (TRUE) {
        # ORF
        if (!is.null(isoData$PTC) &
            'ORF' %in% annotationToAnalyze) {
            # ORF
            isoData$ORF <- 'With ORF'
            isoData$ORF[which(is.na(isoData$PTC))] <- 'Without ORF'
        }
        # PTC
        if (!is.null(isoData$PTC) &
            'NMD_status' %in% annotationToAnalyze) {
            # PTC
            isoData$PTC <-
                ifelse(test = isoData$PTC,
                       yes = 'NMD sensitive',
                       no = 'NMD insensitive')
        } else {
            isoData$PTC <- NULL
        }
        # coding potential
        if (!is.null(isoData$codingPotential)) {
            isoData$codingPotential <-
                ifelse(test = isoData$codingPotential,
                       yes = 'Isoform is coding',
                       no = 'Isoform is non-coding')
        }
        # signal peptide
        if (!is.null(isoData$signal_peptide_identified)) {
            isoData$signal_peptide_identified <-
                ifelse(
                    test = isoData$signal_peptide_identified == 'yes',
                    yes = 'With signal peptide',
                    no = 'Without signal peptide'
                )
        }
        # protein domains
        if (!is.null(isoData$domain_identified)) {
            isoData$domain_identified <-
                ifelse(
                    test = isoData$domain_identified == 'yes',
                    yes = 'With protein domain',
                    no = 'Without protein domain'
                )
        }
        # intron retention
        if (!is.null(isoData$IR)) {
            isoData$IR <-
                ifelse(test = isoData$IR > 0,
                       yes = 'With intron retention',
                       no = 'Without intron retention')
        }
        # switch consequences
        if (!is.null(isoData$switchConsequences)) {
            isoData$switchConsequences <-
                ifelse(
                    test = isoData$switchConsequences,
                    yes = 'With switch consequence',
                    no = 'Without switch consequence'
                )
        }
        # class code
        if (!is.null(isoData$class_code)) {
            isoData$class_code <-
                paste('class code: \"', isoData$class_code, '\"', sep = '')
        }
    }

    ### Prepare for plotting
    if (TRUE) {
        # melt categories
        isoDataMelt <-
            melt(isoData,
                 id.vars = c(
                     'iso_ref', 'isoform_id', 'comparison', 'IF1', 'IF2'
                 ))
        isoDataMelt <-
            isoDataMelt[which(!is.na(isoDataMelt$value)), ]

        # melt IF
        colnames(isoDataMelt)[match(
            c('variable', 'value'), colnames(isoDataMelt)
        )] <- c('category', 'isoform_feature')
        isoDataMelt <-
            melt(
                isoDataMelt,
                id.vars = c(
                    'iso_ref',
                    'isoform_id',
                    'comparison',
                    'category',
                    'isoform_feature'
                )
            )

        # massage category
        isoDataMelt$category <- as.character(isoDataMelt$category)

        isoDataMelt$category <- gsub(
            'PTC',
            'NMD Status',
            isoDataMelt$category
        )
        isoDataMelt$category <-
            gsub('codingPotential',
                 'Coding Potential',
                 isoDataMelt$category)
        isoDataMelt$category <-
            gsub('signal_peptide_identified',
                 'Signal Peptide',
                 isoDataMelt$category)
        isoDataMelt$category <-
            gsub('domain_identified',
                 'Protein Domains',
                 isoDataMelt$category)
        isoDataMelt$category <-
            gsub('IR','Intron Retention',isoDataMelt$category)
        isoDataMelt$category <-
            gsub('switchConsequences',
                 'Switch Consequence',
                 isoDataMelt$category)
        isoDataMelt$category <-
            gsub('class_code','Isoform Class',isoDataMelt$category)

        isoDataMelt$comparison2 <-
            paste(isoDataMelt$comparison, '\n(IF1 vs IF2)', sep = '')
    }

    ### Calculate statistics
    if (TRUE) {
        mySigTest <-
            ddply(
                isoDataMelt,
                .variables = c('comparison', 'category', 'isoform_feature'),
                .drop = TRUE,
                .fun = function(aDF) {
                    data1 <- aDF$value[which(aDF$variable == 'IF1')]
                    data2 <- aDF$value[which(aDF$variable == 'IF2')]

                    myTest <- suppressWarnings(wilcox.test(data1,
                                                           data2))

                    myResult <- data.frame(
                        n = length(data1),
                        medianIF1 = median(data1),
                        medianIF2 = median(data2)
                    )
                    myResult$medianDIF <-
                        myResult$medianIF2 - myResult$medianIF1

                    myResult$wilcoxPval <- myTest$p.value
                    myResult$ymax <- max(c(data1, data2))

                    return(myResult)
                }
            )
        mySigTest$ymax <- mySigTest$ymax * 1.05
        mySigTest$wilcoxQval <-
            p.adjust(mySigTest$wilcoxPval, method = 'fdr')
        mySigTest$significance <-
            sapply(mySigTest$wilcoxQval, function(x)
                evalSig(x, alphas))


        mySigTest <-
            ddply(
                mySigTest,
                .variables = 'category',
                .fun = function(aDF) {
                    aDF$isoform_feature <- factor(aDF$isoform_feature)
                    aDF$idNr <- as.numeric(aDF$isoform_feature)
                    return(aDF)
                }
            )

        mySigTest$comparison2 <-
            paste(mySigTest$comparison, '\n(IF1 vs IF2)', sep = '')

    }

    ### Plot result
    if (plot) {
        # start plot
        if (violinPlot) {
            p1 <- ggplot() +
                geom_violin(
                    data = isoDataMelt,
                    aes(
                        x = isoform_feature,
                        y = value,
                        fill = variable
                    ),
                    scale = 'area'
                ) +
                stat_summary(
                    data = isoDataMelt,
                    aes(
                        x = isoform_feature,
                        y = value,
                        fill = variable
                    ),
                    fun.y = medianQuartile,
                    geom = 'point',
                    position = position_dodge(width = 0.9),
                    size = 2
                )
        } else {
            p1 <- ggplot() +
                geom_boxplot(data = isoDataMelt, aes(
                    x = isoform_feature,
                    y = value,
                    fill = variable
                ))
        }

        # add significance
        p1 <- p1 +
            geom_text(
                data = mySigTest,
                aes(x = isoform_feature, y = ymax, label = significance),
                vjust = -0.2,
                size = localTheme$text$size * 0.3
            ) +
            geom_segment(data = mySigTest, aes(
                x = idNr - 0.25,
                xend = idNr + 0.25,
                y = ymax,
                yend = ymax
            ))

        # build rest of plot
        p1 <- p1 +
            facet_grid(comparison2 ~ category,
                       scales = 'free_x',
                       space = 'free_x') +
            localTheme +
            theme(strip.text.y = element_text(angle = 0)) +
            theme(axis.text.x = element_text(
                angle = -45,
                hjust = 0,
                vjust = 1
            )) +
            scale_fill_discrete(name = NULL) + theme(legend.position = "top") +
            labs(x = 'Isoform feature', y = 'Isoform Usage (IF)') +
            #coord_cartesian(ylim=c(0,1.25))
            coord_cartesian(ylim = c(0, 1.1 + max(c(
                0, 0.02 * (length(unique(
                    mySigTest$comparison2
                )) - 2)
            ))))

        print(p1)
    }

    ### Return result
    if (returnResult) {
        mySigTest2 <-
            mySigTest[, c(
                'comparison',
                'category',
                'isoform_feature',
                'n',
                'medianIF1',
                'medianIF2',
                'medianDIF',
                'wilcoxPval',
                'wilcoxQval',
                'significance'
            )]
        mySigTest2 <-
            mySigTest2[order(mySigTest2$comparison,
                             mySigTest2$category,
                             mySigTest2$isoform_feature), ]
        return(mySigTest2)
    }
}
