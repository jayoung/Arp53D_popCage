library(dplyr)
library(tidyr)

rm(list=ls())

setwd("/Volumes/malik_h/user/jayoung/forOtherPeople/forCourtney/Arp53D_population_cage_experiments/Arp53D_popCage")
source("Arp53D_populationCage_functions.R")


##### read in Courtney's ACTUAL data 
## generations 1 and 5 are not REALLY data, they are dictated by how the experiment was set up

arp53d <- read.table(text="bottle	1	5	10	15	20
1	25	16.69658887	43.814433	63.7119114	63.4011091
2	25	49.8826291	37.704918	66.5427509	65.2747253
3	25	29.13907285	52.6666667	NA	72.2532588", header=TRUE)
colnames(arp53d) <- gsub("X", "gen", colnames(arp53d))
arp53d[,"gen2"] <- 0
arp53d <- arp53d %>% 
    pivot_longer(cols=-bottle, names_to="generation", values_to="percentWThom") %>% 
    mutate(generation=as.integer(gsub("gen","",generation))) %>% 
    mutate(freqWThom=percentWThom/100) %>% 
    arrange(bottle,generation) %>% 
    filter(!is.na(percentWThom))

#### show individual WThom freqs at each generation that represents REAL data (rather than being determined by experimental setup). These numbers are what I'll use to determine the fit of each model
realWTfreqEachGenToFit <- arp53d %>% 
    filter(!generation %in% c(1,2)) %>% 
    select(generation,freqWThom)

realWTfreqEachGenToFit %>% arrange(generation)
# A tibble: 11 x 2
#   generation freqWThom
# 1          5     0.167
# 2          5     0.499
# 3          5     0.291
# 4         10     0.438
# 5         10     0.377
## etc

##### we will only look at the IDEAL case of infinite population size and Mendelian inheritance.  

## set up with our initial genotype freqs (they are imbalanced between males and females, which constrains the genotype frequencies in the first 2 generations)
genotypeFreqsInitial <- list(female=list(WThom=0, het=0, KOhom=1 ), 
                             male=list(WThom=0.5, het=0, KOhom=0.5) )



### now, set up a range of selection regimes, including neutral, that I will use to plug into my simulation. These all assume that the fitness deficiency in a het is half of that of a KOhom. Now a version that is much more fine-grained, with 0.001 intervals between coefficients:
temp <- seq(0,0.6,by=0.001) 
selectionCoefficientsTableFineGrain <- list()
selectionCoefficientsTableFineGrain[[1]] <- data.frame(WThom=rep(1,length(temp)),
                                              het=rep(1,length(temp)),
                                              KOhom=1-temp,
                                              regime="het_equalsWThom")
selectionCoefficientsTableFineGrain[[2]] <- data.frame(WThom=rep(1,length(temp)),
                                              het=1-(temp/2),
                                              KOhom=1-temp,
                                              regime="het_intermediate")
selectionCoefficientsTableFineGrain[[3]] <- data.frame(WThom=rep(1,length(temp)),
                                              het=1-temp,
                                              KOhom=1-temp,
                                              regime="het_equalsKOhom")
selectionCoefficientsTableFineGrain <- do.call("rbind",selectionCoefficientsTableFineGrain)
rm(temp)



######## run all those selective regimes on my experimental setup for only 30 generations, fine-grained selection coefficients:
infinitePopulation_multipleSelectiveRegimesFineGrain <- runModellingOnMultipleSelectiveRegimes(
    selectionTable=selectionCoefficientsTableFineGrain, 
    initialGenotypeFreqs=genotypeFreqsInitial, 
    numGenerations=35)

infinitePopulation_multipleSelectiveRegimesFineGrain_fits <- calculateModelFitsMultipleModels(
    infinitePopulation_multipleSelectiveRegimesFineGrain,
    realDataTable=realWTfreqEachGenToFit)

# show best model use MAE criterion
infinitePopulation_multipleSelectiveRegimesFineGrain_fits %>% 
    filter(fitType=="meanAbsErr") %>% 
    group_by(selectionRegime) %>% 
    filter(myFits == min(myFits))





########  What selective pressure for each regime best matches the real data?

###### get best models for each combination of selection regime and fitTYpe
allSelectionRegimes <- c("het_equalsWThom", "het_intermediate", "het_equalsKOhom")
bestModelsAllRegimesSeveralFitTypes <- lapply( allSelectionRegimes, 
    function(thisSelectionRegime) {
        #allFitTypes <- c("gen15", "gen20", "meanAbsErr")
        allFitTypes <- c("meanAbsErr")
        bestModelsThisRegime <- lapply( allFitTypes, function(thisFitType) {
            #cat("## getting best model for regime",thisSelectionRegime,"and fitType",thisFitType,"\n")
            thisBestModel <- getBestModelResults(infinitePopulation_multipleSelectiveRegimesFineGrain, 
                                infinitePopulation_multipleSelectiveRegimesFineGrain_fits,
                                modelTypeToChoose=thisSelectionRegime, fitTypeToChoose=thisFitType)
            return(thisBestModel)
        } )
        names(bestModelsThisRegime) <- allFitTypes
        return(bestModelsThisRegime)
    })
names(bestModelsAllRegimesSeveralFitTypes) <- allSelectionRegimes


bestModelsAllRegimesSeveralFitTypesCombined <- lapply(
    bestModelsAllRegimesSeveralFitTypes, 
    function(bestModelsOneRegimeSeveralFitTypes) {
        bestModelsOneRegimeSeveralFitTypes <- lapply(names(bestModelsOneRegimeSeveralFitTypes), function(oneFitType) {
            thisBestModel <- bestModelsOneRegimeSeveralFitTypes[[oneFitType]]
            thisBestModel[,"fitType"] <- oneFitType
            return(thisBestModel)
        })
        bestModelsOneRegimeSeveralFitTypes <- do.call("rbind",bestModelsOneRegimeSeveralFitTypes)
        return(bestModelsOneRegimeSeveralFitTypes)
    })
bestModelsAllRegimesSeveralFitTypesCombined <- do.call("rbind",bestModelsAllRegimesSeveralFitTypesCombined)
rownames(bestModelsAllRegimesSeveralFitTypesCombined) <- NULL


### this shows me I DO get the same result as the website, although I name the generations differently: the website starts at generation 0 and I start at generation 1
#infinitePopulation_multipleSelectiveRegimes_website %>% 
#    filter(generation==21 & selectionKOhom==0.84 & selectionRegime=="het_intermediate")

save.image(file="Rdata_files/Arp53D_populationCage_modellingOutput.Rdata")


