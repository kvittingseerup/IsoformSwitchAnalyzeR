### TO DO ###


###################################################
### Introduce alternative splicing from spliceR ###
###################################################

### Done ###
#  -  Implement analysis of differences in consequence gain/loss (like AS)
#       - incl documentation
#  - include analyzeAlternativSplicing() as the main AS analysis function.
#       - replace the current analyzeIntronRetention() with this analyzeIntronRetention() which is just a warapper for analyzeAlternativSplicing()
#       - Include all funcions needed from spliceR - done
#  - remove spliceR dependency
#  - move and rename extractGenomeWideAnalysis() rename to extractGenomeWideConsequenceAnalysis()
#  - R Documentation for the new functions
#  - Check build sucess
#  - AS names must be changed in annotation - not functions
#  - analyzeAlternativSplicing() adds the entry "AlternativeSplicingAnalysis" to the switch list
#        analyzeSwitchConsequences() needs to look under both "AlternativeSplicingAnalysis" and "IntronRetentionAnalysis" for the IR analysis
#  - Update vignette
#       - High level workflow
#           - Include the summary splicing
#       - Detailed workflow
#           - Update of splicing section
#           - Update of [Genome Wide Analysis]
# - Update example data
#  - Add examples to R documentation of new functions
# - Improve splicing summary - they are slow - done
#  - The getCDS function is problematc? - done
# Add functionality for comparing different comparisons
#   - Compare enrichments via fisher.test
#   - sepearte section in other workflows wich summarizes the functions needed
#   - analyzeAlternativeSplicing instead of analyzeAlternativSplicing()
#   - Color of significant / nonsignificant in enrichment plots
#   - UpsetR plots
# - SummarizeSplicing
#       Why are combinedfactions not summing to 1.
#       Why are MME factions 1.
# - Better x-axis text in enrichment
# Update switchTest to only perform correction if all conditions agree



##### OTHER

### Add mean gene/isoform expression + mean IF and allow filtering on those

### Update vignette so importing data from Kallisto/Salmon is the default is the one shown in quick overview

# Update isoformToGeneExp works more like tximport
#   Because counts should be calculated from abundances at gene level

# Add a run_info part to the swichList which allows for easy backgracking
#   (and wanrs you when you are using different cutoffs)

