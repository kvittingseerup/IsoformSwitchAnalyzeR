### Helper function for p-value calibration
pvcEstimateSigma2 <- function(
    ### The implementation of the p-value correction found in Ferguson et al P-value calibration for multiple testing problems in genomics ###
    # This is a modified version of the R function published with the article. The modification was simply to split the function into two functions.
    # This modified version only contains the parts for esimating and returning the Sigma2.
    # R code: http://www.math.bas.bg/~palejev/pvc/pvc.r
    # vignette : http://www.math.bas.bg/~palejev/pvc/

    ### Arguments
    pvalues,        # the original p-values
    pvaluesUse=NA,  # indicates which p-values would be used when estimating the likelihood. Possible values: either NA (default) or a boolean vector of the same length as pvalues. If NA, all p-values would be used when estimating the likelihood.If boolean vector of the same length as pvalues, only the elements of pvalues sliced by pvalueUse would be used when estimating the likelihood
    startsigma2=1,  # starting value of sigma2 in conditional likelihood and EM algorithm.
    condlik = TRUE, # set to FALSE if EM-algorithm based calibration is desired
    # conditional likelihood
    neff=NA,        # a prior estimate of the number of effective genes. Only useful when plotting the conditional likelihood
    plotl=FALSE,    # A logic indicting if a plot of the likelihood is desired, when running the conditional likelihood method
    le=0.1,          # lower bound for inclusion of p-values into the likelihood when running conditional likelihood method
    ue=1,           # upper bound for inclusion of p-values into the likelihood when running conditional likelihood method
    # EM implementation
    starteps=0.01,  # starting value for epsilon parameter in EM algorithm
    startmu=1,      # starting value for mu parameter in EM algorithm
    tol=.00001      # convergence criterion for EM algorithm

    ### Output
    # sigma2 - estimate of sigma2
) {
    Np <- length(pvalues)
    if((length(pvaluesUse)==1) && is.na(pvaluesUse)) {
        pvaluesUse <- rep(TRUE, Np)}
    if(length(pvaluesUse)!=Np) return("pvaluesUse not specified correctly")
    myindexes <- c(1:length(pvalues))[!is.na(pvalues)]
    pvalues <- pvalues[myindexes]
    pvaluesUse <- pvaluesUse[myindexes]
    prob_cal_needed <- NA

    ### Conditional likelyhood
    if(condlik){
        indexes <- c(1:length(pvalues))[pvalues < ue & pvalues >= le & pvaluesUse]
        tv <- qnorm(1-pvalues[indexes]/2)
        N <- length(tv)
        c1 <- qnorm(1-le/2)
        c2 <- qnorm(1-ue/2)
        simple_l <- function(sigma2){
            -(N/2)*log(sigma2)+sum(-1*tv^2/(2*sigma2))-
                N*log(pnorm(c1/sqrt(sigma2))-pnorm(c2/sqrt(sigma2)))}
        if(!is.na(neff)) simple_l <- function(sigma2){-(neff/2)*log(sigma2)+
                (neff/N)*sum(-1*tv^2/(2*sigma2))-neff*
                log(pnorm(c1/sqrt(sigma2))-pnorm(c2/sqrt(sigma2)))}
        sigma2new <- optim(par=startsigma2,method="L-BFGS-B",function(x){
            -1*simple_l(x)},lower=0.01, upper=10)$par
        if(plotl) {
            sigma_possible <- seq(from = 0.2, to = 2, by = 0.0001)
            like_vals <- sapply(sigma_possible, function(x){simple_l(x)})
            like_vals <- like_vals - max(like_vals)
            cum_area <- cumsum(exp(like_vals))/sum(exp(like_vals))
            prob_cal_needed <- max(cum_area[sigma_possible < 0.9])
            plot(sigma_possible, exp(like_vals)/(.0001*sum(exp(like_vals))),
                 type="l",main="", xlab=expression(sigma^2),ylab = "")
            title(expression(paste("Normalized likelihood for ",sigma^2)))
            S <- .0001*sum(exp(like_vals))
            for(j in 1:length(cum_area[sigma_possible < 0.9])) lines(
                x=rep(sigma_possible[j],2),y=c(0,exp(like_vals[j])/S),
                col="palevioletred")
            text(x=.75,y=.5,labels=paste(
                "area = ",round(prob_cal_needed,2),sep=""))
        }
    }

    ### EM algorithm
    if(!condlik){
        indexes <- c(1:length(pvalues))[pvalues<1 & pvaluesUse]
        tv <- qnorm(1-pvalues[indexes]/2)
        tv[tv=="Inf"]=max(tv[tv!="Inf"])
        N <- length(tv)

        sigma2old <- startsigma2
        currenterror <- 1
        muold <- startmu
        epsold <- starteps
        while(currenterror > tol){
            ddiff <- dnorm((tv-muold)/sqrt(sigma2old))+dnorm((-1*tv-muold)/sqrt(sigma2old))
            dnormal <- 2*dnorm(tv/sqrt(sigma2old))
            probs <- epsold*ddiff/((epsold*ddiff)+(1-epsold)*dnormal)
            a1 <- epsold*dnorm((tv-muold)/sqrt(sigma2old))/((epsold*ddiff)+(1-epsold)*dnormal)
            a2 <- epsold*dnorm((-1*tv-muold)/sqrt(sigma2old))/((epsold*ddiff)+(1-epsold)*dnormal)
            munew <- weighted.mean(c(tv,-1*tv), w=c(a1,a2))
            sigma2new <- (1/N)*sum((1-probs)*tv^2 + a1*(tv-munew)^2+a2*(-1*tv-munew)^2)
            epsnew <- mean(probs)
            currenterror <- max(abs(epsnew-epsold),abs(sigma2new-sigma2old),abs(munew-muold))
            sigma2old <- sigma2new
            muold <- munew
            epsold <- epsnew

        }
        postprobs= rep(NA,length(pvalues))
        postprobs[indexes]=probs
        if(length(myindexes)<Np){
            postprobs1 <- rep(NA,Np)
            postprobs1[myindexes] <- postprobs
            postprobs <- postprobs1
        }
        logL <- -1*sum((1-probs)*tv^2)/(2*sigma2new)-1*sum(a1*(tv-munew)^2)/(2*sigma2new)-1*sum(a2*(-1*tv-munew)^2)/(2*sigma2new)-(N/2)*log(2*pi*sigma2new)+log(epsnew)*sum(probs)+log(1-epsnew)*(N-sum(probs))
    }

    return(sigma2new)

}

pvcApplySigma <- function(
    ### The implementation of the p-value correction found in Ferguson et al P-value calibration for multiple testing problems in genomics ###
    # This is a modified version of the R function published with the article. The modification was simply to split the function into two functions.
    # This modified version only contains the parts for using the estimated sigma2 to callibrated the p-values.
    # R code: http://www.math.bas.bg/~palejev/pvc/pvc.r
    # vignette : http://www.math.bas.bg/~palejev/pvc/

    ### Arguments
    pvalues,        # the original p-values
    sigma2new

    ### Output
    # calibrated p-values
) {
    Np <- length(pvalues)
    myindexes <- c(1:length(pvalues))[!is.na(pvalues)]
    pvalues <- pvalues[myindexes]

    ### transform the original pvalues.
    indexes <- c(1:length(pvalues))[pvalues<1]
    tv <- qnorm(1-pvalues[indexes]/2)
    tv[tv=="Inf"]=max(tv[tv!="Inf"])
    N <- length(tv)
    pvalues= rep(1,length(pvalues))
    pvalues[indexes] <- pchisq(tv^2/sigma2new,0,lower.tail=FALSE,df=1)
    if(length(myindexes)<Np){
        pvalues1 <- rep(NA,Np)
        pvalues1[myindexes] <- pvalues
        pvalues <- pvalues1
    }

    return(pvalues)
}


### Test for isoform switching
isoformSwitchTest <- function(
    switchAnalyzeRlist,
    alpha = 0.05,
    dIFcutoff = 0.1,
    reduceToSwitchingGenes = TRUE,
    calibratePvalues = TRUE,
    showProgress = FALSE,
    quiet = FALSE
) {
    ### Test input
    if (TRUE) {
        ### Check cufflinks version (if data originates from cufflinks)
        if (grepl('^cufflinks_', switchAnalyzeRlist$sourceId)) {
            checkVersionFail <- function(versionVector,
                                         minVersionVector) {
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
                as.integer(strsplit(
                    gsub('cufflinks_', '', switchAnalyzeRlist$sourceId),
                    '\\.'
                )[[1]])
            if (checkVersionFail(cuffVersionDeconstructed, c(2, 2, 1))) {
                stop(
                    paste(
                        'The version of cufflinks you have used is outdated.',
                        'An error in the estimations of expression standard',
                        'errors was not corrected untill',
                        'cufflinks v2.2.1. Since the detection of isoform',
                        'switches using this function',
                        'relies on this standard error estimat, isoform',
                        'switches cannot be tested from this',
                        'versions of cufflinks. Please upgrade cufflinks',
                        'to >=version 2.2.1 or newer and try again.',
                        sep = ' '
                    )
                )
            }
        }

        ### Tjek number of replicates
        if (any(switchAnalyzeRlist$conditions$nrReplicates == 1)) {
            stop(paste(
                'This function relies on biological replicates to estimate',
                'the untertainties in the data and can therefore',
                'not be used in comparisons without replicates'
            ))
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

        ### Tjek arguments
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist')        {
            stop(paste(
                'The object supplied to \'switchAnalyzeRlist\'',
                'must be a \'switchAnalyzeRlist\''
            ))
        }
        if (!is.logical(reduceToSwitchingGenes))  {
            stop(paste(
                'The argument supplied to \'reduceToSwitchingGenes\'',
                'must be an a logic'
            ))
        }
        if (dIFcutoff < 0 | dIFcutoff > 1) {
            stop('The dIFcutoff must be in the interval [0,1].')
        }

    }

    if (showProgress &  !quiet) {
        progressBar <- 'text'
    } else {
        progressBar <- 'none'
    }

    ### Extract and filter expression data
    if (!quiet) {
        message('Step 1 of 3: Filtering for eligible data...')
    }
    if (TRUE) {
        ### Extract data
        columnsToExtract <-
            c(
                'iso_ref',
                'gene_ref',
                "isoform_id",
                "gene_id",
                "condition_1",
                "condition_2",
                "gene_value_1",
                "gene_value_2",
                'gene_stderr_1',
                'gene_stderr_2',
                "iso_value_1",
                "iso_value_2",
                'iso_stderr_1',
                'iso_stderr_2',
                'IF1',
                'IF2',
                'dIF'
            )
        myDiffData <-
            switchAnalyzeRlist$isoformFeatures[, columnsToExtract]

        ### Filter for null variance
        myDiffData <-
            myDiffData[which(myDiffData$iso_stderr_1 != 0 &
                                 myDiffData$iso_stderr_2 != 0),]
        myDiffData <-
            myDiffData[which(!is.na(myDiffData$gene_stderr_1) |
                                 !is.na(myDiffData$gene_stderr_2)), ]
        if (!nrow(myDiffData)) {
            stop(paste(
                'No isoforms were eligible for isoform switch',
                'analysis after removing zero-variance genes'
            ))
        }

        ### Extract genes with multiple isoforms
        # get genes with multiple isoforms
        geneIsoOverview <-
            unique(myDiffData[, c('isoform_id', 'gene_id')])
        isoList <-
            split(geneIsoOverview$isoform_id, f = geneIsoOverview$gene_id)
        genesOfInterest <-
            names(isoList[which(sapply(isoList, length) > 1)])
        # do the filtering
        myDiffData <-
            myDiffData[which(myDiffData$gene_id %in% genesOfInterest), ]
        if (!nrow(myDiffData)) {
            stop(paste(
                'No isoforms were eligible for isoform switch',
                'analysis after removing single-isoform genes'
            ))
        }

        ### Add replcate number
        myDiffData$nrReplicates_1 <-
            switchAnalyzeRlist$conditions$nrReplicates[match(
                myDiffData$condition_1 ,
                switchAnalyzeRlist$conditions$condition)]
        myDiffData$nrReplicates_2 <-
            switchAnalyzeRlist$conditions$nrReplicates[match(
                myDiffData$condition_2 ,
                switchAnalyzeRlist$conditions$condition)]

        myDiffData$gene_cv1 <-
            (myDiffData$gene_stderr_1 * sqrt(myDiffData$nrReplicates_1)) /
            myDiffData$gene_value_1
        myDiffData$gene_cv2 <-
            (myDiffData$gene_stderr_2 * sqrt(myDiffData$nrReplicates_2)) /
            myDiffData$gene_value_2

        ### Calculate CI
        confidenceInterval <-
            0.95 # hardcoded since this is not subject to change
        myDiffData$gene_lower_CI_1 <-
            myDiffData$gene_value_1 - myDiffData$gene_stderr_1 *
            qt(confidenceInterval / 2 + .5,
               df = myDiffData$nrReplicates_1 - 1)
        myDiffData$gene_lower_CI_2 <- myDiffData$gene_value_2 -
            myDiffData$gene_stderr_2 *
            qt(confidenceInterval /2 + .5, df = myDiffData$nrReplicates_2 - 1)

        ### Filter on CV and CI - this is nesseary
        # for the implemnted version of Fieller's theorem
        myDiffData2 <-
            myDiffData[which(myDiffData$gene_lower_CI_1 > 0 &
                                 myDiffData$gene_lower_CI_2 > 0), ] # This will also remove the genes with gene NA in standard error.
        myDiffData2 <- myDiffData2[which(
            myDiffData2$gene_stderr_1 * sqrt(myDiffData2$nrReplicates_1) <
                myDiffData2$gene_value_1 / 2 &
                myDiffData2$gene_stderr_2 * sqrt(myDiffData2$nrReplicates_2) <
                myDiffData2$gene_value_2 / 2
        ), ]

        if (!nrow(myDiffData2)) {
            stop(paste(
                'No isoforms were eligible for isoform switch analysis',
                'after removing genes expressed to close to zero'
            ))
        }
    }

    ### Calculate isoform fraction variance
    if (TRUE) {
        ### Calculate expression variannce
        myDiffData2$gene_var_1 <-
            (myDiffData2$gene_stderr_1 * sqrt(myDiffData2$nrReplicates_1)) ^ 2
        myDiffData2$gene_var_2 <-
            (myDiffData2$gene_stderr_2 * sqrt(myDiffData2$nrReplicates_2)) ^ 2

        myDiffData2$iso_var_1 <-
            (myDiffData2$iso_stderr_1 * sqrt(myDiffData2$nrReplicates_1)) ^ 2
        myDiffData2$iso_var_2 <-
            (myDiffData2$iso_stderr_2 * sqrt(myDiffData2$nrReplicates_2)) ^ 2

        ### calculate var of isoform fraction
        myDiffData2$IF_var_1 <- (1 / myDiffData2$gene_value_1 ^ 2) *
            (
                myDiffData2$iso_var_1 +
                    (myDiffData2$IF1 ^ 2 * myDiffData2$gene_var_1) -
                    (2 * myDiffData2$IF1 * myDiffData2$iso_var_1)
            )

        myDiffData2$IF_var_2 <- (1 / myDiffData2$gene_value_2 ^ 2) *
            (
                myDiffData2$iso_var_2 +
                    (myDiffData2$IF2 ^ 2 * myDiffData2$gene_var_2) -
                    (2 * myDiffData2$IF2 * myDiffData2$iso_var_2)
            )


        myDiffData2 <-
            myDiffData2[which(myDiffData2$IF_var_1 > 0 &
                                  myDiffData2$IF_var_2 > 0), ]

        ### Print message
        localFraction <-
            round(
                nrow(myDiffData2) /
                    nrow(switchAnalyzeRlist$isoformFeatures) * 100,
                  digits = 2)
        if (!quiet) {
            message(
                paste(
                    'Found',
                    nrow(myDiffData2),
                    '(',
                    localFraction ,
                    '%)',
                    'of isoform comparisons eligable for switch analysis',
                    sep = ' '
                )
            )
        }
    }

    ### Do statistical test of IF differences
    if (!quiet) {
        message('Step 2 of 3: Testing isoform usage...')
    }
    if (TRUE) {
        ### Split in comparisons (since different
        # comparisons might have different numbert of replicates)
        myDiffData2List <-
            split(myDiffData2, f = myDiffData2[, c(
                'condition_1', 'condition_2')
                ], drop = TRUE)

        ### Apply over each condition, do test, callibrate and correct p-values
        myDiffData3 <-
            do.call(
                rbind,
                llply(
                    .data = myDiffData2List,
                    .progress = progressBar,
                    .fun = function(aDF) {
                        # aDF <- myDiffData2List[[1]]
                        ### Calculate standard error of dIF
                        aDF$dIF_std_err <-
                            sqrt(
                                (aDF$IF_var_1 / aDF$nrReplicates_1) +
                                    (aDF$IF_var_2 / aDF$nrReplicates_2)
                            )

                        ### Calculate test statistics
                        aDF$t_statistics <- aDF$dIF / aDF$dIF_std_err

                        ### Calulate degrees of freedome
                        aDF$deg_free <-
                            (
                                (aDF$IF_var_1 / aDF$nrReplicates_1) +
                                    (aDF$IF_var_2 / aDF$nrReplicates_2)
                            ) ^ 2    /    (((aDF$IF_var_1 ^ 2) / (
                                aDF$nrReplicates_1 ^ 2 *
                                    (aDF$nrReplicates_1 - 1)
                            ))  +  ((aDF$IF_var_2 ^ 2) / (
                                aDF$nrReplicates_2 ^ 2 *
                                    (aDF$nrReplicates_2 - 1)
                            )))

                        ### Calculate p value
                        aDF$p_value <-
                            2 * pt(abs(aDF$t_statistics),
                                   df = aDF$deg_free,
                                   lower.tail = FALSE) # The 2* to make it two tailed

                        aDF <- aDF[which(!is.na(aDF$p_value)),]
                        if (nrow(aDF) == 0) {
                            return(NULL)
                        }

                        ### Perform p-value callibration
                        # Create subset of highly expressed isoforms
                        isoExpCutoff <-
                            max(c(1, quantile(
                                c(aDF$iso_value_1, aDF$iso_value_2) ,
                                probs = 0.50
                            )))

                        highlyExpressedIsoforms <- aDF[which(
                            aDF$iso_value_1  > isoExpCutoff  &
                            aDF$iso_value_2  > isoExpCutoff
                            ),]

                        # estimate sigma squared from highly expressed data
                        estimatedSigma2 <-
                            pvcEstimateSigma2(highlyExpressedIsoforms$p_value)

                        # perform p-value callibration and adjusment
                        if (calibratePvalues &
                            estimatedSigma2 < 0.9) {
                            # in accordance with article advice
                            aDF$calibrated_p_values <-
                                pvcApplySigma(pvalues = aDF$p_value,
                                              sigma2new = estimatedSigma2)
                            aDF$isoform_switch_q_value <-
                                p.adjust(aDF$calibrated_p_values, method = 'BH')
                        } else {
                            aDF$calibrated_p_values <- NA
                            aDF$isoform_switch_q_value <-
                                p.adjust(aDF$p_value , method = 'BH')
                        }

                        return(aDF)
                    }
                )
            )

        ### Extrapolate to gene level
        geneQlevel <- sapply(
            X = split(
                myDiffData3$isoform_switch_q_value,
                f = myDiffData3$gene_ref
            ),
            FUN = function(x) {
                min(c(1, x), na.rm = TRUE)
            }
        )
        myDiffData3$gene_switch_q_value <-
            geneQlevel[match(myDiffData3$gene_ref , names(geneQlevel))]
        rownames(myDiffData3) <- NULL
    }

    ### Report results
    if (!quiet) {
        message('Step 3 of 3: Preparing output...')
    }
    if (TRUE) {
        ### Overwrite previous results
        switchAnalyzeRlist$isoformFeatures$gene_switch_q_value <- NA
        switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value <-
            NA

        ### Add significance to isoform_feature tab
        switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value <-
            myDiffData3$isoform_switch_q_value[match(
                switchAnalyzeRlist$isoformFeatures$iso_ref,
                myDiffData3$iso_ref
            )]
        switchAnalyzeRlist$isoformFeatures$gene_switch_q_value  <-
            myDiffData3$gene_switch_q_value[match(
                switchAnalyzeRlist$isoformFeatures$gene_ref,
                myDiffData3$gene_ref
            )]

        ### Reduce data to significnat if needed
        if (reduceToSwitchingGenes) {
            isoResTest <-
                any(!is.na(
                    switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
                ))
            if (isoResTest) {
                combinedGeneIDsToKeep <- myDiffData3$gene_ref[which(
                    myDiffData3$isoform_switch_q_value < alpha &
                        abs(myDiffData3$dIF) > dIFcutoff
                )]
            } else {
                combinedGeneIDsToKeep <- myDiffData3$gene_ref[which(
                    myDiffData3$gene_switch_q_value    < alpha &
                        abs(myDiffData3$dIF) > dIFcutoff
                )]
            }
            if (length(combinedGeneIDsToKeep) == 0) {
                stop(paste(
                    'No signifcant switches were found with the supplied',
                    'cutoffs whereby we cannot reduce the',
                    'switchAnalyzeRlist to only significant genes'
                ))
            }

            myDiffData3 <-
                myDiffData3[which(myDiffData3$gene_ref %in%
                                      combinedGeneIDsToKeep), ]
            switchAnalyzeRlist <-
                subsetSwitchAnalyzeRlist(
                    switchAnalyzeRlist,
                    switchAnalyzeRlist$isoformFeatures$gene_ref %in%
                        combinedGeneIDsToKeep
                )
        }

        ### Add detailed switch analysis to switchAnalyzeRlist
        columnsToExtract <-
            c(
                'iso_ref',
                'gene_ref',
                "isoform_id",
                "gene_id",
                "condition_1",
                "condition_2",
                'nrReplicates_1',
                'nrReplicates_2',
                'IF1',
                'IF2',
                'IF_var_1',
                'IF_var_2',
                'dIF',
                'dIF_std_err',
                't_statistics',
                'deg_free',
                'p_value',
                'calibrated_p_values',
                'isoform_switch_q_value',
                'gene_switch_q_value'
            )

        myDiffData4 <- myDiffData3[, columnsToExtract]
        rownames(myDiffData4) <- NULL

        switchAnalyzeRlist$isoformSwitchAnalysis <- myDiffData4

        ### test callibration status
        if (calibratePvalues) {
            if (length(unique(
                extractCalibrationStatus(switchAnalyzeRlist)$pvalsAreCalibrated
            )) != 1) {
                warning(
                    paste(
                        'Note that not all comparison was callibrated',
                        '- which means',
                        'they cannot be directly be compared as they have',
                        'been treated differently.',
                        'See ?extractCalibrationStatus for more information or',
                        'run isoformSwitchTest() with calibratePvalues=FALSE.',
                        sep = ' '
                    )
                )
            }
        }

    }

    return(switchAnalyzeRlist)
}


### Test via DRIMSeq
isoformSwitchTestDRIMSeq <- function(
    switchAnalyzeRlist,
    alpha = 0.05,
    dIFcutoff = 0.1,
    testIntegration = 'isoform_only',
    reduceToSwitchingGenes = TRUE,
    #nCores=1,
    dmPrecisionArgs = list(),
    dmFitArgs = list(),
    dmTestArgs = list(),
    showProgress = TRUE,
    quiet = FALSE
) {
    ### Test data
    if (TRUE) {
        ### Tjek arguments
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist')        {
            stop(paste(
                'The object supplied to \'switchAnalyzeRlist\'',
                'must be a \'switchAnalyzeRlist\''
            ))
        }
        if (!is.logical(reduceToSwitchingGenes))  {
            stop(paste(
                'The argument supplied to \'reduceToSwitchingGenes\'',
                'must be an a logic'
            ))
        }
        if (dIFcutoff < 0 | dIFcutoff > 1) {
            stop('The dIFcutoff must be in the interval [0,1].')
        }
        if (reduceToSwitchingGenes) {
            if (alpha < 0 |
                alpha > 1) {
                stop('The alpha parameter must be between 0 and 1 ([0,1]).')
            }
            if (alpha > 0.05) {
                warning(paste(
                    'Most journals and scientists consider an alpha',
                    'larger than 0.05 untrustworthy. We therefore recommend',
                    'using alpha values smaller than or queal to 0.05'
                ))
            }
        }
        if (is.null(switchAnalyzeRlist$isoformCountMatrix)) {
            stop(paste(
                'An isoform replicate count matrix is nesseary for using',
                'the DRIMSeq approach - please reinitalize the',
                'swithcAnalyzeRlist object with one of the import*()',
                'functions and try again.'
            ))
        }

        if (length(testIntegration) > 1) {
            stop('The \'testIntegration\' argument must be of length 1')
        }
        if (!testIntegration %in% c('isoform_only', 'gene_only', 'intersect')) {
            stop(paste(
                'The \'testIntegration\' argument must be one',
                'of \'isoform_only\', \'gene_only\', \'intersect\''
            ))
        }

    }

    ### Extract comparison
    comaprisonsToMake <- unique(switchAnalyzeRlist$isoformFeatures[, c(
        'condition_1', 'condition_2'
    )])

    if (showProgress &  !quiet &  nrow(comaprisonsToMake) > 1) {
        progressBar <- 'text'
    } else {
        progressBar <- 'none'
    }

    ### For each comparison do the test
    # Note : when contrasts are implemented the dmFit() can be moved out and only calculated once
    if (!quiet) {
        message('Step 1 of 3: Analyzing comparisons (this may take a while)...')
    }
    resultOfPairwiseTest <- myListToDf(
        dlply(
            .data = comaprisonsToMake,
            .variables = c('condition_1', 'condition_2'),
            .progress = progressBar,
            .fun = function(aDF) {
                # aDF <- comaprisonsToMake[1,]

                ### Makde design matrix
                if (TRUE) {
                    ### Extract local data
                    localDesign <-
                        switchAnalyzeRlist$designMatrix[which(
                            switchAnalyzeRlist$designMatrix$condition %in%
                                c(aDF$condition_1, aDF$condition_2)
                        ), ]
                    rownames(localDesign) <- localDesign$sampleID

                    ### Ensure correct intercept
                    localDesign$condition <- factor(as.character(
                        localDesign$condition),
                        levels = c(c(
                            aDF$condition_1, aDF$condition_2
                        )))

                    ### Make model matrix
                    # formula
                    localFormula <- paste('~ ', colnames(localDesign)[2])
                    if (ncol(localDesign) > 2) {
                        localFormula <- paste(
                            localFormula,
                            '+',
                            paste(
                                colnames(localDesign)[3:ncol(localDesign)],
                                collapse = '+'
                            )
                        )
                    }
                    localFormula <- as.formula(localFormula)

                    # model
                    localModel <-
                        model.matrix(localFormula, data = localDesign)
                    colnames(localModel)[2] <- 'coefOfInterest'

                    ### Massage design
                    localDesign2 <- localDesign[, 1:2]
                    colnames(localDesign2) <- c('sample_id', 'group')
                }

                ### Extract corresponding count matrix
                if (TRUE) {
                    localCM <- switchAnalyzeRlist$isoformCountMatrix[, c(
                        'isoform_id', rownames(localModel))]

                    ### Add gene_ref for easy backward tracking
                    localCM$gene_id <-
                        switchAnalyzeRlist$isoformFeatures$gene_ref[match(
                            paste0(
                                localCM$isoform_id,
                                aDF$condition_1,
                                aDF$condition_2
                            ),
                            paste0(
                                switchAnalyzeRlist$isoformFeatures$isoform_id,
                                switchAnalyzeRlist$isoformFeatures$condition_1,
                                switchAnalyzeRlist$isoformFeatures$condition_2
                            )
                        )]
                    colnames(localCM)[1] <- 'feature_id'

                    ### Subset to genes expressed
                    localCM <- localCM[which(!is.na(localCM$gene_id)), ]

                    ### Massage CM
                    localCM <- localCM[,c(
                        'feature_id', 'gene_id',
                        setdiff(colnames(localCM), c('feature_id', 'gene_id'))
                    )]
                }

                ### Do test
                if (TRUE) {
                    ### Construct DM object
                    suppressMessages(localDm <-dmDSdata(
                        counts = localCM, samples = localDesign2
                    ))

                    ### Calculate precision
                    # add data arguments to argument list
                    dmPrecisionArgs$x      <- localDm
                    dmPrecisionArgs$design <- localModel

                    # use argument list to run function
                    suppressMessages(localDmPrec <- do.call(
                        what = dmPrecision, args = dmPrecisionArgs
                    ))

                    ### Calculate fit
                    # add data arguments to argument list
                    dmFitArgs$x        <- localDmPrec
                    dmFitArgs$design   <- localModel
                    dmFitArgs$bb_model <- TRUE

                    # use argument list to run function
                    suppressMessages(localDmFit <- do.call(
                        what = dmFit, args = dmFitArgs
                    ))


                    ### Make test
                    # add data arguments to argument list
                    dmTestArgs$x    <- localDmFit
                    dmTestArgs$coef <- 'coefOfInterest'

                    # use argument list to run function
                    suppressMessages(localDmTest <- do.call(
                        what = dmTest,
                        args = dmTestArgs
                    ))

                    ### Extract result
                    localRes <- merge(
                        x = results(localDmTest, level = "feature")[, c(
                            'feature_id',
                            'gene_id',
                            'lr',
                            'df',
                            'pvalue',
                            'adj_pvalue'
                        )],
                        y = results(localDmTest)[, c(
                            'gene_id', 'lr', 'df', 'pvalue', 'adj_pvalue'
                            )],
                        by = 'gene_id',
                        suffixes = c(".iso", ".gene")
                    )
                }

                return(localRes)
            }
        )
    )

    ### Massage result
    if (!quiet) {
        message('Step 2 of 3: Preparing output...')
    }
    if (TRUE) {
        ### Remove NAs
        resultOfPairwiseTest <-
            resultOfPairwiseTest[which(
                !is.na(resultOfPairwiseTest$adj_pvalue.iso)
                ), ]

        isoInd <-
            which(grepl('iso$', colnames(resultOfPairwiseTest)))
        geneInd <-
            which(grepl('gene$', colnames(resultOfPairwiseTest)))
        colnames(resultOfPairwiseTest) <-
            gsub('\\.gene|\\.iso', '', colnames(resultOfPairwiseTest))

        colnames(resultOfPairwiseTest)[isoInd] <- paste0(
            'iso_', colnames(resultOfPairwiseTest)[isoInd])
        colnames(resultOfPairwiseTest)[geneInd] <- paste0(
            'gene_', colnames(resultOfPairwiseTest)[geneInd])
        colnames(resultOfPairwiseTest) <-
            gsub('adj_pvalue',
                 'q_value',
                 colnames(resultOfPairwiseTest))
        colnames(resultOfPairwiseTest) <-
            gsub('pvalue', 'p_value', colnames(resultOfPairwiseTest))

        colnames(resultOfPairwiseTest) <-
            gsub('feature_id',
                 'isoform_id',
                 colnames(resultOfPairwiseTest))


        myDiff <-
            setdiff(colnames(resultOfPairwiseTest),
                    c('isoform_id', 'gene_id'))
        resultOfPairwiseTest <- resultOfPairwiseTest[, c(
            'isoform_id', 'gene_id',
            myDiff[which(grepl('^gene', myDiff))],
            myDiff[which(grepl('^iso', myDiff))]
        )]
        colnames(resultOfPairwiseTest)[2] <- 'gene_ref'

        resultOfPairwiseTest$isoform_id <-
            switchAnalyzeRlist$isoformFeatures$iso_ref[match(
                paste0(
                    resultOfPairwiseTest$isoform_id,
                    resultOfPairwiseTest$gene_ref
                ),
                paste0(
                    switchAnalyzeRlist$isoformFeatures$isoform_id,
                    switchAnalyzeRlist$isoformFeatures$gene_ref
                )
            )]
        colnames(resultOfPairwiseTest)[1] <- 'iso_ref'
    }

    ### Add result to switchAnalyzeRlist
    if (!quiet) {
        message('Step 3 of 3: Adding result to switchAnalyzeRlist...')
    }
    if (TRUE) {
        ### Overwrite previous results
        switchAnalyzeRlist$isoformFeatures$gene_switch_q_value <- NA
        switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value <-
            NA

        ## Interpret the p-values via the testIntegration argument
        if (testIntegration == 'isoform_only') {
            ### summarize to gene level
            geneQlevel <- sapply(
                X = split(
                    resultOfPairwiseTest$iso_q_value,
                    f = resultOfPairwiseTest$gene_ref
                ),
                FUN = function(x) {
                    min(c(1, x), na.rm = TRUE)
                }
            )
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value <-
                geneQlevel[match(switchAnalyzeRlist$isoformFeatures$gene_ref,
                                 names(geneQlevel))]

            ### Isoform level
            switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value <-
                resultOfPairwiseTest$iso_q_value[match(
                    switchAnalyzeRlist$isoformFeatures$iso_ref,
                    resultOfPairwiseTest$iso_ref
                )]

        } else if (testIntegration == 'gene_only') {
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value <-
                resultOfPairwiseTest$gene_q_value[match(
                    switchAnalyzeRlist$isoformFeatures$iso_ref,
                    resultOfPairwiseTest$iso_ref
                )]

            switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value <-
                NA

        } else if (testIntegration == 'intersect') {
            localRes <-
                resultOfPairwiseTest[which(
                    resultOfPairwiseTest$iso_q_value >=
                        resultOfPairwiseTest$gene_p_value
                ), ]

            ### summarize to gene level
            geneQlevel <- sapply(
                X = split(localRes$iso_q_value, f = localRes$gene_ref),
                FUN = function(x)
                    min(c(1, x), na.rm = TRUE)
            )
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value <-
                geneQlevel[match(switchAnalyzeRlist$isoformFeatures$gene_ref,
                                 names(geneQlevel))]

            ### Isoform level
            switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value <-
                localRes$iso_q_value[match(
                    switchAnalyzeRlist$isoformFeatures$iso_ref,
                    localRes$iso_ref
                )]

        } else {
            stop(paste(
                'Somthing went wrong with the integraion of p-values',
                '- please contact developer with a reproducible example'
            ))
        }

        ### Add the full analysis
        switchAnalyzeRlist$isoformSwitchAnalysis <-
            resultOfPairwiseTest

    }

    ### Print status
    if (!quiet) {
        myN <-
            length(unique(switchAnalyzeRlist$isoformSwitchAnalysis$gene_ref))
        myFrac <-
            myN / length(unique(
                switchAnalyzeRlist$isoformFeatures$gene_ref
                )) * 100
        message(
            paste(
                'An isoform switch analysis was performed for ',
                myN,
                ' genes (',
                round(myFrac, digits = 1),
                '%).',
                sep = ''
            )
        )
    }

    ### Reduce to genes with switches
    if (reduceToSwitchingGenes) {
        isoResTest <-
            any(!is.na(
                switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
            ))
        if (isoResTest) {
            combinedGeneIDsToKeep <-
                switchAnalyzeRlist$isoformFeatures$gene_ref[which(
                    switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value <
                        alpha &
                        abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
                )]
        } else {
            combinedGeneIDsToKeep <-
                switchAnalyzeRlist$isoformFeatures$gene_ref[which(
                    switchAnalyzeRlist$isoformFeatures$gene_switch_q_value <
                        alpha &
                        abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
                )]
        }
        if (length(combinedGeneIDsToKeep) == 0) {
            stop(paste(
                'No signifcant switches were found with the supplied cutoffs',
                'whereby we cannot reduce the switchAnalyzeRlist to only',
                'significant genes'
            ))
        }

        switchAnalyzeRlist <-
            subsetSwitchAnalyzeRlist(
                switchAnalyzeRlist,
                switchAnalyzeRlist$isoformFeatures$gene_ref %in% combinedGeneIDsToKeep
            )
    }


    if (!quiet) {
        message('Done')
    }
    return(switchAnalyzeRlist)
}



### Summarize switching
extractSwitchSummary <- function(
    switchAnalyzeRlist,
    filterForConsequences = FALSE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    includeCombined = nrow(unique(
        switchAnalyzeRlist$isoformFeatures[, c('condition_1', 'condition_1')]
    )) > 1
) {
    ### Test input
    if (TRUE) {
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(paste(
                'The object supplied to \'switchAnalyzeRlist\' must',
                'be a \'switchAnalyzeRlist\''
            ))
        }
        if (!any(!is.na(
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value
        ))) {
            stop(paste(
                'The analsis of isoform switching must be performed before',
                'functional consequences can be analyzed. Please run ?isoformSwitchTest and try again.'
            ))
        }
        if (filterForConsequences) {
            if (!'switchConsequencesGene' %in%
                colnames(switchAnalyzeRlist$isoformFeatures)) {
                stop(paste(
                    'The switchAnalyzeRlist does not contain isoform',
                    'switching analysis. Please run the',
                    '\'isoformSwitchTest\' function first.'
                ))
            }
        }
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
        if (dIFcutoff < 0 | dIFcutoff > 1) {
            stop('The dIFcutoff must be in the interval [0,1].')
        }
    }

    backUpDf <-
        unique(switchAnalyzeRlist$isoformFeatures[, c(
            'condition_1', 'condition_2'
        )])
    backUpDf <-
        data.frame(
            Comparison = paste(
                backUpDf$condition_1,
                backUpDf$condition_2, sep = ' vs '),
            nrIsoforms = 0,
            nrGenes = 0,
            stringsAsFactors = FALSE
        )

    ### Extract data needed
    columnsToExtract <-
        c(
            'isoform_id',
            'gene_id',
            'condition_1',
            'condition_2',
            'dIF',
            'isoform_switch_q_value',
            'gene_switch_q_value',
            'switchConsequencesGene'
        )

    isoResTest <-
        any(!is.na(
            switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
        ))
    if (isoResTest) {
        dataDF <- switchAnalyzeRlist$isoformFeatures[which(
            switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value < alpha &
                abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
        ),
        na.omit(match(
            columnsToExtract,
            colnames(switchAnalyzeRlist$isoformFeatures)
        ))]
    } else {
        dataDF <- switchAnalyzeRlist$isoformFeatures[which(
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value    < alpha &
                abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
        ),
        na.omit(match(
            columnsToExtract,
            colnames(switchAnalyzeRlist$isoformFeatures)
        ))]
    }
    if (nrow(dataDF) == 0) {
        return(backUpDf)
    }

    if (filterForConsequences) {
        dataDF <- dataDF[which(dataDF$switchConsequencesGene), ]
        if (nrow(dataDF) == 0) {
            return(backUpDf)
        }
    }

    ### Summarize pr comparison
    dataList <-
        split(dataDF,
              f = paste(dataDF$condition_1, dataDF$condition_2, sep = ' vs '))

    if (length(dataList) > 1 | includeCombined) {
        dataList$combined <- dataDF
    }

    isoResTest <-
        any(!is.na(
            switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
        ))
    myNumbers <- ldply(dataList, function(aDF) {
        if (isoResTest) {
            sigGenes <-
                unique(aDF$gene_id   [which(aDF$isoform_switch_q_value < alpha &
                                                abs(aDF$dIF) > dIFcutoff)])
            sigIso   <-
                unique(aDF$isoform_id[which(aDF$isoform_switch_q_value < alpha &
                                                abs(aDF$dIF) > dIFcutoff)])
            return(data.frame(
                nrIsoforms = length(sigIso),
                nrGenes = length(sigGenes)
            ))
        } else {
            sigGenes <-
                unique(aDF$gene_id   [which(aDF$gene_switch_q_value < alpha &
                                                abs(aDF$dIF) > dIFcutoff)])
            sigIso   <-
                unique(aDF$isoform_id[which(aDF$gene_switch_q_value < alpha &
                                                abs(aDF$dIF) > dIFcutoff)])
            return(data.frame(
                nrIsoforms = length(sigIso),
                nrGenes = length(sigGenes)
            ))
        }


    })
    colnames(myNumbers)[1] <- 'Comparison'
    return(myNumbers)
}

### Extract
extractTopSwitches <- function(
    switchAnalyzeRlist,
    filterForConsequences = FALSE,
    extractGenes = TRUE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    n = 10,
    inEachComparison = FALSE,
    sortByQvals = TRUE
) {
    ### Test input
    if (TRUE) {
        if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
            stop(paste(
                'The object supplied to \'switchAnalyzeRlist\'',
                'must be a \'switchAnalyzeRlist\''
            ))
        }
        if (!any(!is.na(
            switchAnalyzeRlist$isoformFeatures$gene_switch_q_value
        ))) {
            stop(paste(
                'The analsis of isoform switching must be performed before',
                'functional consequences can be analyzed.',
                'Please run ?isoformSwitchTest and try again.'
            ))
        }
        if (filterForConsequences) {
            if (!'switchConsequencesGene' %in%
                colnames(switchAnalyzeRlist$isoformFeatures)) {
                stop(paste(
                    'The switchAnalyzeRlist does not contain isoform',
                    'switching analysis. Please run the \'isoformSwitchTest\'',
                    'function first.'
                ))
            }
        }
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
        if (dIFcutoff < 0 | dIFcutoff > 1) {
            stop('The dIFcutoff must be in the interval [0,1].')
        }
        if (!is.logical(extractGenes)) {
            stop('The extractGenes argument supplied must be a logical')
        }
        if (!is.logical(inEachComparison)) {
            stop('The inEachComparison argument supplied must be a logical')
        }
        if (!is.logical(sortByQvals)) {
            stop('The sortByQvals argument supplied must be a logical')
        }
    }

    if (extractGenes) {
        columnsToExtract <-
            c(
                'gene_id',
                'gene_name',
                'condition_1',
                'condition_2',
                'gene_switch_q_value',
                'switchConsequencesGene',
                'dIF'
            )
        isoResTest <-
            any(!is.na(
                switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
            ))
        if (isoResTest) {
            dataDF <- unique(switchAnalyzeRlist$isoformFeatures[which(
                switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value <
                    alpha &
                    abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
            ),
            na.omit(match(
                columnsToExtract,
                colnames(switchAnalyzeRlist$isoformFeatures)
            ))])
        } else {
            dataDF <- unique(switchAnalyzeRlist$isoformFeatures[which(
                switchAnalyzeRlist$isoformFeatures$gene_switch_q_value <
                    alpha &
                    abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
            ),
            na.omit(match(
                columnsToExtract,
                colnames(switchAnalyzeRlist$isoformFeatures)
            ))])
        }
        if (nrow(dataDF) == 0) {
            stop('No significant switching genes were found.')
        }

        if (filterForConsequences) {
            dataDF <- dataDF[which(dataDF$switchConsequencesGene), ]
        }
        if (nrow(dataDF) == 0) {
            stop('No significant switching genes consequences were found.')
        }

        ### Sort data
        dataDF$comparison <-
            paste(dataDF$condition_1, '_vs_', dataDF$condition_2, sep = '')
        if (is.na(sortByQvals)) {
            dataDF2 <- dataDF
            rownames(dataDF2) <- NULL
        } else if (sortByQvals) {
            dataDF2 <-
                dataDF[sort.list(
                    dataDF$gene_switch_q_value,
                    decreasing = FALSE
                ), ]
            rownames(dataDF2) <- NULL
        } else {
            dataDF$combinedID <- paste(dataDF$gene_id,
                                       dataDF$condition_1,
                                       dataDF$condition_2,
                                       sep = '_')

            ### Calculate combined dIF value
            combinedDif <-
                split(abs(dataDF$dIF), f = dataDF$combinedID)
            combinedDif <- sapply(combinedDif, sum)

            ### Add to df
            dataDF$combinedDIF <-
                combinedDif[match(dataDF$combinedID , names(combinedDif))]
            dataDF$combinedID <- NULL

            dataDF2 <-
                dataDF[sort.list(dataDF$combinedDIF, decreasing = TRUE), ]
            rownames(dataDF2) <- NULL
        }

        ### reduce to collumns wanted
        columnsToExtract <-
            c(
                'gene_id',
                'gene_name',
                'condition_1',
                'condition_2',
                'gene_switch_q_value',
                'combinedDIF',
                'switchConsequencesGene'
            )
        dataDF2 <-
            unique(dataDF2[, na.omit(
                match(columnsToExtract, colnames(dataDF))
            )])


        ### Reduce to the number wanted
        if (!is.na(n)) {
            if (inEachComparison) {
                dataDF2$comparison <-
                    paste(dataDF2$condition_1,
                          '_vs_',
                          dataDF2$condition_2,
                          sep = '')
            } else {
                dataDF2$comparison <- 'AllCombined'
            }

            dataDF2 <-
                ddply(
                    .data = dataDF2,
                    .variables = 'comparison',
                    .fun = function(aDF) {
                        if (nrow(dataDF2) < n) {
                            if (filterForConsequences) {
                                warning(paste(
                                    'Less than n genes with significant',
                                    'switches and consequences were found.',
                                    'Returning those.'
                                ))
                            } else {
                                warning(paste(
                                    'Less than n genes genes with significant',
                                    'switches were found. Returning those.'
                                ))
                            }
                            n2 <- nrow(dataDF)
                        } else {
                            n2 <- n
                        }

                        return(aDF[1:n2, ])
                    }
                )

            dataDF2$comparison <- NULL
        }

        return(dataDF2)
    }

    ### Extract data to retun
    if (!extractGenes) {
        columnsToExtract <-
            c(
                'isoform_id',
                'gene_id',
                'gene_name',
                'condition_1',
                'condition_2',
                'iso_significant',
                'IF1',
                'IF2',
                'dIF',
                'isoform_switch_q_value',
                'switchConsequencesGene'
            )
        isoResTest <-
            any(!is.na(
                switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
            ))
        if (isoResTest) {
            dataDF <- unique(switchAnalyzeRlist$isoformFeatures[which(
                switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value <
                    alpha &
                    abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
            ),
            na.omit(match(
                columnsToExtract,
                colnames(switchAnalyzeRlist$isoformFeatures)
            ))])
        } else {
            dataDF <- unique(switchAnalyzeRlist$isoformFeatures[which(
                switchAnalyzeRlist$isoformFeatures$gene_switch_q_value <
                    alpha &
                    abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
            ),
            na.omit(match(
                columnsToExtract,
                colnames(switchAnalyzeRlist$isoformFeatures)
            ))])
        }
        if (nrow(dataDF) == 0) {
            stop('No significant switching isoforms were found.')
        }

        if (filterForConsequences) {
            dataDF <- dataDF[which(dataDF$switchConsequencesGene), ]
        }
        if (nrow(dataDF) == 0) {
            stop(
                'No significant switching isoforms with consequences were found.'
            )
        }

        ### Sort data
        dataDF$comparison <-
            paste(dataDF$condition_1, '_vs_', dataDF$condition_2, sep = '')
        if (is.na(sortByQvals)) {
            dataDF2 <- dataDF
            rownames(dataDF2) <- NULL
        } else if (sortByQvals) {
            dataDF2 <-
                dataDF[sort.list(
                    dataDF$isoform_switch_q_value,
                    decreasing = FALSE
                ), ]
            rownames(dataDF2) <- NULL
        } else {
            dataDF2 <- dataDF[sort.list(abs(dataDF$dIF), decreasing = TRUE), ]
            rownames(dataDF2) <- NULL
        }

        ### Reduce to the number wanted
        if (!is.na(n)) {
            if (inEachComparison) {
                dataDF2$comparison <-
                    paste(dataDF2$condition_1,
                          '_vs_',
                          dataDF2$condition_2,
                          sep = '')
            } else {
                dataDF2$comparison <- 'AllCombined'
            }

            dataDF2 <-
                ddply(
                    .data = dataDF2,
                    .variables = 'comparison',
                    .fun = function(aDF) {
                        if (nrow(dataDF2) < n) {
                            if (filterForConsequences) {
                                warning(paste(
                                    'Less than n genes with significant',
                                    'switches and consequences were found.',
                                    'Returning those.'
                                ))
                            } else {
                                warning(paste(
                                    'Less than n genes genes with significant',
                                    'switches were found. Returning those.'
                                ))
                            }
                            n2 <- nrow(dataDF)
                        } else {
                            n2 <- n
                        }

                        return(aDF[1:n2, ])
                    }
                )

            dataDF2$comparison <- NULL
        }

        dataDF2$IF1 <- round(dataDF2$IF1, digits = 3)
        dataDF2$IF2 <- round(dataDF2$IF2, digits = 3)
        dataDF2$dIF <- round(dataDF2$dIF, digits = 3)

        return(dataDF2)
    }

}


extractCalibrationStatus <- function(
    switchAnalyzeRlist
) {
    ### Test input
    if (TRUE) {
        if (is.null(switchAnalyzeRlist$isoformSwitchAnalysis)) {
            stop(paste(
                'This function can only analyse switchAnalyzeRlist which have',
                'been analyzed with the \'isoformSwitchTest\'',
                'function. Please run \'isoformSwitchTest()\' and try again.'
            ))
        }
    }

    ### Extract status
    callibratedPvals <- split(
        switchAnalyzeRlist$isoformSwitchAnalysis$calibrated_p_values,
        f = paste(
            switchAnalyzeRlist$isoformSwitchAnalysis$condition_1,
            ' vs ',
            switchAnalyzeRlist$isoformSwitchAnalysis$condition_2,
            sep = ''
        )
    )

    myStatus <- ldply(
        callibratedPvals,
        .fun = function(aVec) {
            data.frame(pvalsAreCalibrated = any(!is.na(aVec)))
        }
    )
    colnames(myStatus)[1] <- 'Comparison'

    return(myStatus)
}
