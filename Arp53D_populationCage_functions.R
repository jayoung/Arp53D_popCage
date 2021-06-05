
########### functions for mating infinite size population

##### mendelTable: set up a table that shows, for every possible pairwise combination of parent genotypes, the offspring genotype frequencies under Mendelian inheritance:
genotypePossibilities <- c("WThom", "het", "KOhom")
mendelTable <- data.frame(femaleID=rep(genotypePossibilities, each=3), 
                          maleID=rep(genotypePossibilities, 3))

nextGenColumnNames <- paste("nextGen_",genotypePossibilities, sep="")
mendelTable[,nextGenColumnNames] <- NA

# mother is WThom
mendelTable[which(mendelTable[,"femaleID"]=="WThom" & mendelTable[,"maleID"]=="WThom"),
            nextGenColumnNames] <- c(1,0,0)
mendelTable[which(mendelTable[,"femaleID"]=="WThom" & mendelTable[,"maleID"]=="het"),
            nextGenColumnNames] <- c(0.5,0.5,0)
mendelTable[which(mendelTable[,"femaleID"]=="WThom" & mendelTable[,"maleID"]=="KOhom"),
            nextGenColumnNames] <- c(0,1,0)

# mother is het
mendelTable[which(mendelTable[,"femaleID"]=="het" & mendelTable[,"maleID"]=="WThom"),
            nextGenColumnNames] <- c(0.5,0.5,0)
mendelTable[which(mendelTable[,"femaleID"]=="het" & mendelTable[,"maleID"]=="het"),
            nextGenColumnNames] <- c(0.25,0.5,0.25)
mendelTable[which(mendelTable[,"femaleID"]=="het" & mendelTable[,"maleID"]=="KOhom"),
            nextGenColumnNames] <- c(0,0.5,0.5)

# mother is KOhom
mendelTable[which(mendelTable[,"femaleID"]=="KOhom" & mendelTable[,"maleID"]=="WThom"),
            nextGenColumnNames] <- c(0,1,0)
mendelTable[which(mendelTable[,"femaleID"]=="KOhom" & mendelTable[,"maleID"]=="het"),
            nextGenColumnNames] <- c(0,0.5,0.5)
mendelTable[which(mendelTable[,"femaleID"]=="KOhom" & mendelTable[,"maleID"]=="KOhom"),
            nextGenColumnNames] <- c(0,0,1)
rm(nextGenColumnNames)


####### mateInfiniteSizePopulation - a function that takes as input 
## (a) parental genotypes (b) table of offspring genotype freqs for each mating pair type
## it first figures out frequencies of each possible pairwise combination of mating pair genotypes (P-Gmat x Gpat)
## then for each mating pair type, it figures out frequency of offspring genotypes (O-wt, O-het, O-ko)
## it takes the product of those frequencies (PxO) and sums across all mating pair types to get the overall frequency of each genotype in the next generation
# (a) the first input, the parental genotypes, is a list of lists: first level of the list is female/male, each of which is a list with three elements whose values are genotype frequencies (female frequencies must sum to 1, same for male)
mateInfiniteSizePopulation <- function(genotypeFreqList, mendel=mendelTable) {
    ## make sure genotypes look OK
    verifyGenotypeFreqList(genotypeFreqList)

    # we will mate each combination of parents separately and add up the resulting offspring
    nextGen <- data.frame(femaleID= rep(names(genotypeFreqList[["female"]]),each=3), 
                          femaleFreq=rep(unlist(genotypeFreqList[["female"]], use.names=FALSE),each=3),
                          maleID= rep(names(genotypeFreqList[["male"]]), 3), 
                          maleFreq=rep(unlist(genotypeFreqList[["male"]], use.names=FALSE), 3))
    nextGen[,"parentCombFreq"] <- nextGen[,"femaleFreq"] * nextGen[,"maleFreq"]
    
    if(!identical (mendel[,"femaleID"], nextGen[,"femaleID"])) {
        cat("\n\nERROR - something is wrong in mateInfiniteSizePopulation - nextGen table and mendel table have different row orders\n\n")
        return(list(mendel=mendel, nextGen=nextGen))
    }
    if(!identical (mendel[,"maleID"], nextGen[,"maleID"])) {
        cat("\n\nERROR - something is wrong in mateInfiniteSizePopulation - nextGen table and mendel table have different row orders\n\n")
        return(list(mendel=mendel, nextGen=nextGen))
    }
    
    nextGenColumnNames <- colnames(mendel)[3:5]
    # check row totals of mendel table add up to 1
    if (sum(apply(mendel[,nextGenColumnNames], 1, sum) != 1) > 0) {
        badRows <- which(apply(mendel[,nextGenColumnNames], 1, sum) != 1) 
        cat("\n\nERROR in mateInfiniteSizePopulation - one or more rows of the mendel table do not sum to 1 - here they are:\n\n")
        print(mendel[badRows,])
        cat("\n\n")
        stop("\n\nERROR in mateInfiniteSizePopulation - see above\n\n")
    }
    
    # get offspring freqs for each row using mendelian genetics
    nextGen[,nextGenColumnNames] <- NA
    for (i in 1:9) {
        nextGen[i,nextGenColumnNames] <- nextGen[i,"parentCombFreq"] * mendel[i,nextGenColumnNames]
    }
    
    genotypeFreqTotals <- apply(nextGen[,nextGenColumnNames], 2, sum)
    names(genotypeFreqTotals) <- gsub("nextGen_","",names(genotypeFreqTotals))
    return(list(table=nextGen, totals=genotypeFreqTotals))
}

#### verifyGenotypeFreqList - a function to check that the genotype frequencies list looks as it should
verifyGenotypeFreqList <- function(genotypeFreqList, anyGenotype=genotypePossibilities, errorTolerance=0.0000001) {
    if (length(genotypeFreqList) != 2) {
        stop("\n\nERROR in verifyGenotypeFreqList - it should be a list object with two elements (female and male)\n\n")
    }
    if (sum(c("female","male") %in% names(genotypeFreqList)) != 2) {
        stop("\n\nERROR in verifyGenotypeFreqList - elements of the genotype frequency list should be named female and male\n\n")
    }
    if (length(genotypeFreqList[["female"]]) != 3) {
        stop("\n\nERROR in verifyGenotypeFreqList - female genotype frequencies should have 3 elements\n\n")
    } 
    if (length(genotypeFreqList[["male"]]) != 3) {
        stop("\n\nERROR in verifyGenotypeFreqList - male genotype frequencies should have 3 elements\n\n")
    } 
    genotypeFreqSums <- lapply(genotypeFreqList, function(x) { sum(unlist(x)) })
    
    ## I think rounding errors mean these don't always QUITE sum to 1, so I'll be a little tolerant
    if (abs(1-genotypeFreqSums[["female"]]) > errorTolerance) {
        stop("\n\nERROR in verifyGenotypeFreqList - female genotype frequencies should sum to 1\n\n")
    } 
    if (abs(1-genotypeFreqSums[["male"]]) > errorTolerance) {
        stop("\n\nERROR in verifyGenotypeFreqList - male genotype frequencies should sum to 1\n\n")
    } 
    if (sum(anyGenotype %in% names(genotypeFreqList[["female"]])) != 3) {
        stop("\n\nERROR in verifyGenotypeFreqList - elements of the female genotype frequencies list should be named",anyGenotype,"\n\n")
    }
    if (sum(anyGenotype %in% names(genotypeFreqList[["male"]])) != 3) {
        stop("\n\nERROR in verifyGenotypeFreqList - elements of the male genotype frequencies list should be named",anyGenotype,"\n\n")
    }
    return(TRUE)
}


#### getAlleleFreqsGivenGenotypeFreqs: a function to get the two allele frequencies, given a simple numeric vector of the three genotype frequencies
getAlleleFreqsGivenGenotypeFreqs <- function(genotypes, expectedGenotypes=genotypePossibilities) {
    if(class(genotypes)=="list") { genotypes <- unlist(genotypes) }
    ## make sure the genotype freqs have the expected names:
    if(!identical(sort(expectedGenotypes), sort(names(genotypes)))) {
        stop("\n\nERROR in getAlleleFreqsGivenGenotypeFreqs function - I do not see all the expected genotype classes in the genotypes you supplied: should have:",expectedGenotypes,"\n\n")
    }
    alleleFreqs <- numeric()
    alleleFreqs["WT"] <- ((2 * genotypes["WThom"])  +  genotypes["het"]) / 2
    alleleFreqs["KO"] <- ((2 * genotypes["KOhom"])  +  genotypes["het"]) / 2
    return(alleleFreqs)
}


#### runInfinitePopulationModelling - a function to iterate mateInfiniteSizePopulation over multiple generations, perhaps applying selection to the offspring genotypes at each generation
runInfinitePopulationModelling <- function(genotypeFreqsInitial, 
                                           numGenerations=NULL, 
                                           selectionCoefficients=NULL,
                                           debug=FALSE) {
    
    initialPopulation <- list()
    initialPopulation[[1]] <- list()
    initialPopulation[[1]][["totals"]] <- genotypeFreqsInitial
    initialPopulation[[1]][["alleleFreqs"]] <- apply(sapply(genotypeFreqsInitial, function(x) { 
        getAlleleFreqsGivenGenotypeFreqs(unlist(x)) 
    }), 1, mean)
    
    ## at first I tried normalizing the selectionCoefficients so that they add up to 1, but I realised it is better to do this on the selected genotypes: that works better especially if some genotypes have frequency=0
    #if(!is.null(selectionCoefficients)) {
    #    normalizedSelectionCoefficients <- selectionCoefficients / mean(selectionCoefficients)
    #}
    outputPopulations <- initialPopulation
    for (i in 2:numGenerations) {
        if (debug) {cat("generation",i,"\n")}
        prevGenGenotypeFreqs <- outputPopulations[[(i-1)]][["totals"]]
        # when previous generation is 1, prevGenGenotypeFreqs looks a bit different, as we separately specified male and female
        if( !"female" %in% names(prevGenGenotypeFreqs) ) {
            prevGenGenotypeFreqs <- list(female=prevGenGenotypeFreqs, male=prevGenGenotypeFreqs) 
        }
        outputPopulations[[i]] <- mateInfiniteSizePopulation(prevGenGenotypeFreqs)
        
        if(!is.null(selectionCoefficients)) {
            ## we use selectionCoefficients to adjust the genotype totals of the offspring.
            if(!identical(names(selectionCoefficients), 
                          names(outputPopulations[[i]][["totals"]]))) {
                stop("\n\nERROR in runInfinitePopulationModelling - the names of the selection coefficients are not what we expected\n\n")
            }
            outputPopulations[[i]][["totalsBeforeSelection"]] <- outputPopulations[[i]][["totals"]]
            
            newGenotypeTotals <- selectionCoefficients * outputPopulations[[i]][["totals"]]
            # normalize so that newGenotypeTotals still adds up to 1
            newGenotypeTotals <- newGenotypeTotals / sum(newGenotypeTotals)
            outputPopulations[[i]][["totals"]] <- newGenotypeTotals
        }
        outputPopulations[[i]][["alleleFreqs"]] <- getAlleleFreqsGivenGenotypeFreqs(outputPopulations[[i]][["totals"]])
        rm(prevGenGenotypeFreqs)
    }
    return(outputPopulations)
}


### summarizeInfinitePopulationModelling - a function to summarize genotype and allele frequencies from runInfinitePopulationModelling
summarizeInfinitePopulationModelling <- function(modellingOutput) {
    summaryTable <- data.frame(generation=1:length(modellingOutput))
    summaryTable <- cbind(summaryTable, 
                          as.data.frame(t(sapply(modellingOutput, "[[", "alleleFreqs"))))
    genotypeFreqs <- as.data.frame(t(sapply(modellingOutput, function(x) {
        g <- x[["totals"]]
        # the unusual case where we have separate male and female freqs (we assume males and females are at equal freq):
        if("female" %in% names(g)) { g <- (unlist(g[["female"]]) + unlist(g[["male"]]))/2 } 
        return(g)
    })))
    summaryTable <- cbind(summaryTable, genotypeFreqs)
    return(summaryTable)
} 


#### runModellingOnMultipleSelectiveRegimes is a function to run runInfinitePopulationModelling and summarizeInfinitePopulationModelling on a number of selective regimes. 
# I only keep genotype and allele freq totals, not the other information that runInfinitePopulationModelling yields at each generation
runModellingOnMultipleSelectiveRegimes <- function(selectionTable=NULL, 
                                                   initialGenotypeFreqs=NULL, 
                                                   numGenerations=NULL, 
                                                   debug=FALSE) {
    if(is.null(selectionTable)) {
        stop("\n\nERROR in runModellingOnMultipleSelectiveRegimes - must supply a table of selective pressures using selectionTable argument\n\n")
    }
    if(is.null(initialGenotypeFreqs)) {
        stop("\n\nERROR in runModellingOnMultipleSelectiveRegimes - must supply initial genotype frequencies using initialGenotypeFreqs argument\n\n")
    }
    if(is.null(numGenerations)) {
        stop("\n\nERROR in runModellingOnMultipleSelectiveRegimes - must supply number of generations using numGenerations argument\n\n")
    }
    resultsEachSelectiveRegime <- lapply(1:dim(selectionTable)[1], function(i){
        thisModel <- runInfinitePopulationModelling(
            initialGenotypeFreqs, 
            numGenerations=numGenerations, 
            selectionCoefficients=unlist(selectionTable[i,c("WThom","het","KOhom")]),
            debug=debug)
        summarizedResults <- summarizeInfinitePopulationModelling(thisModel)
        summarizedResults[,"selectionWThom"] <- selectionTable[i,"WThom"]
        summarizedResults[,"selectionHet"] <- selectionTable[i,"het"]
        summarizedResults[,"selectionKOhom"] <- selectionTable[i,"KOhom"]
        summarizedResults[,"selectionRegime"] <- selectionTable[i,"regime"]
        return(summarizedResults)
    })
    resultsEachSelectiveRegime <- do.call("rbind",resultsEachSelectiveRegime)
    resultsEachSelectiveRegime[,"selectionRegime"] <- factor(resultsEachSelectiveRegime[,"selectionRegime"], 
                                                             levels=c("het_equalsWThom", "het_intermediate", "het_equalsKOhom"))
    return(resultsEachSelectiveRegime)
}


######## a couple of functions to check the fit to the real data for each model. 

##### calculateModelFitsOneModel gets the MAE score of a single model
## I use MAE (mean absolute error) to understand how well a model fits the real data, considering all data points supplied in realDataTable. Small numbers of MAE are good!
## I was initially playing around with fitting to the mean at every generation, or the mean at a single generation, but I realized MAE considering all data points is probably better
## we only consider the fit of WThom frequencies, because in our real experiment that's the data we actually have. 
calculateModelFitsOneModel <- function(oneModelResults,realDataTable) {
    allFits <- list()
    ### calculate overall fits using mean absolute error of model versus realDataTable
    individualDiffs <- numeric()
    for (thisGen in unique(realDataTable$generation)) {
        thisGenModelResults <- oneModelResults %>% 
            filter(generation==thisGen)
        if(dim(thisGenModelResults)[1] > 1) {
            stop("\n\nERROR in calculateFitEachModel - there is more than one row in the model results from generation",thisGen,"\n\n")
        }  
        if(dim(thisGenModelResults)[1] < 1) {
            stop("\n\nERROR in calculateFitEachModel - there are no rows in the model results from generation",thisGen,"\n\n")
        }  
        thisGenModelResults <- thisGenModelResults[1,"WThom"]
        thisGenRealResults <- realDataTable %>% filter(generation==thisGen) %>% select(freqWThom) %>% unlist()
        thisGenModelResults <- rep(thisGenModelResults, length(thisGenRealResults))
        individualDiffs <- c(individualDiffs, abs(thisGenModelResults - thisGenRealResults) )
    }
    allFits[["meanAbsErr"]] <- mean(individualDiffs) 
    allFits <- unlist(allFits)
    return(allFits)
}

###### calculateModelFitsMultipleModels is a function to split modelling results table into results from individual models, then get the fit for each model to the real data
calculateModelFitsMultipleModels <- function(allModelsResults, 
                                             realDataTable) {
    # first, split the results table into results for each model separately.  
    # Together, selectionRegime and selectionKOhom should specify a unique selection regime. If they don't, I'll get an error from calculateModelFitsOneModel because there'll be >1 row for each generation
    eachModelList <- split( allModelsResults, list(allModelsResults$selectionRegime, allModelsResults$selectionKOhom))
    # for each modelling setup, get the fit
    eachModelFits <- lapply(eachModelList, function(x) {
        myFitsDF <- x[1,c("selectionWThom", "selectionHet", "selectionKOhom", "selectionRegime")]
        myFits <- calculateModelFitsOneModel(x, realDataTable=realDataTable)
        myFitsDF[,"fitType"] <- names(myFits)
        myFitsDF <- cbind(myFitsDF, myFits)
        return(myFitsDF)
    })    
    eachModelFits <- do.call("rbind",eachModelFits)
    row.names(eachModelFits) <- NULL
    return(eachModelFits)
}


##### given a bunch of modelling results and their first, choose the best one for a given selectionRegime (e.g. het_equalsWThom) and fitType, which is likely meanAbsErr, but in previous versions I was fitting to the mean at generation 15 ('gen15') etc
getBestModelResults <- function(allModelResults, 
                                allModelResultsFitTable, 
                                modelTypeToChoose, 
                                fitTypeToChoose="meanAbsErr") {
    bestModel <- allModelResultsFitTable %>% 
        filter(fitType==fitTypeToChoose & selectionRegime==modelTypeToChoose) %>% 
        filter(myFits == min(myFits))
    # check there's a single best model
    if(dim(bestModel)[1]>1) {
        stop("\n\nERROR in getModelResultsForSingleSelectionRegime - found >1 best model - don't know how to proceed\n\n")
    }
    if(dim(bestModel)[1]<1) {
        stop("\n\nERROR in getModelResultsForSingleSelectionRegime - found NO best model - don't know how to proceed\n\n")
    }
    # now get results for that best model
    bestModelResults <- allModelResults %>% 
        filter(selectionRegime == modelTypeToChoose) %>% 
        filter(selectionWThom == as.numeric(bestModel$selectionWThom)) %>% 
        filter(selectionHet == as.numeric(bestModel$selectionHet)) %>% 
        filter(selectionKOhom == as.numeric(bestModel$selectionKOhom))
    
    # add the fitType and the fitScore as I may want it later:
    bestModelResults[,"fitType"] <- fitTypeToChoose
    bestModelResults[,"fitScore"] <- bestModel[1,"myFits"]
    return(bestModelResults)
}

