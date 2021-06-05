library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(scales)
library(grid)
library(gridExtra)

rm(list=ls())

setwd("/Volumes/malik_h/user/jayoung/forOtherPeople/forCourtney/Arp53D_population_cage_experiments/Arp53D_popCage")

load("Rdata_files/Arp53D_populationCage_modellingOutput.Rdata")
source("Arp53D_populationCage_functions.R")


####### make a plot for the paper:

#### get text labels that show the parameters of the best models
getModelCoefficientsAsLabel_v2 <- function(modelResults) {
    paste("best model:",
          "\n  Fwt = ",modelResults[1,"selectionWThom"],
          "\n Fhet = ",modelResults[1,"selectionHet"],
          "\n  Fko = ",modelResults[1,"selectionKOhom"], 
          "\n  MAE = ",round(modelResults[1,"fitScore"],digits=3), sep="")
}
bestModelLabels_allRegimesFitTypes_v2 <- lapply(bestModelsAllRegimesSeveralFitTypes, function(x) {
    lapply(x, getModelCoefficientsAsLabel_v2)
})
bestModelLabels_allRegimesFitTypes_v2_df <- lapply(names(bestModelLabels_allRegimesFitTypes_v2), function(thisSelectionRegime) {
    theseLabels <- bestModelLabels_allRegimesFitTypes_v2[[thisSelectionRegime]]
    theseLabels_df <- data.frame(fitType=names(theseLabels), label=unlist(theseLabels,use.names = FALSE))
    theseLabels_df[,"selectionRegime"] <- thisSelectionRegime
    return(theseLabels_df)
})
bestModelLabels_allRegimesFitTypes_v2_df <- do.call("rbind",bestModelLabels_allRegimesFitTypes_v2_df) %>% 
    select(selectionRegime,fitType,label) %>% # reorder columns
    mutate(generation=0,WThom=0.85) %>% # add columns that will specify label position on plots
    mutate(selectionRegime=factor(selectionRegime, levels=allSelectionRegimes))


facetNameReplacements <- c(
    "het_equalsWThom"  = "Fhet = Fwt",
    "het_intermediate" = "Fhet = intermediate",
    "het_equalsKOhom"  = "Fhet = Fko"
)





plotFitnessModelling_allRegimes_justMAE <- infinitePopulation_multipleSelectiveRegimesFineGrain %>% 
    filter( ( (selectionKOhom*1000) %% 25)==0 ) %>%  ## this is so I don't plot every single increment of 0.001 for selectionKOhom (instead, I plot increments of 0.025).## I thought I could use modulo to filter for increments, but there is something weird about floating point arithmetic that means it doesn't work.  See https://stackoverflow.com/questions/13614749/modulus-bug-in-r
    # e.g. 1 %% 0.2  should be 0, but on my Mac R thinks it is 0.2
    # so instead I will multiply by 1000 first, and then take the modulo
    ggplot(aes(x=generation, y=WThom)) +
    geom_line(aes(group=selectionKOhom, colour=selectionKOhom)) +
    facet_grid(cols=vars(selectionRegime), labeller=labeller(selectionRegime=facetNameReplacements)) +
    theme_classic() +
    coord_cartesian(xlim=c(0,30)) + 
    scale_colour_distiller(palette = "Spectral", direction=1, 
                           guide = guide_colourbar(title="Fko")) +
    labs(x="Generation", y="Freq WT homozygotes") +
    ## the gray line showing best model
    geom_line(data=(bestModelsAllRegimesSeveralFitTypesCombined %>% filter(fitType=="meanAbsErr")), 
              aes(x=generation,y=WThom), color="gray",lwd=2) +
    
    ## the REAL data:
    geom_point(data=arp53d, aes(x=generation,y=freqWThom)) +
    ## old - I was including means at each generation
    #stat_summary(data=arp53d, aes(x=generation,y=freqWThom), fun="mean", 
    #             geom="point", shape=5, size=3, 
    #             colour = "black", fill="white") +
    ## dashed lines for the REAL data
    geom_line(data=arp53d, aes(x=generation,y=freqWThom, group=bottle), lty=2) +
    ## gray boxes showing coefficients for best models
    geom_label(data=(bestModelLabels_allRegimesFitTypes_v2_df %>% filter(fitType=="meanAbsErr")), 
               aes(label=gsub("fitScore","MAE",label)), 
               colour="gray60", hjust = 0, #vjust = 1, label.padding=unit(0.25, "lines"),
               size=1.75, family = "mono") + 
    theme(panel.spacing = unit(2, "lines"),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 12)
    )


plotFitnessModelling_allRegimes_justMAE


ggsave(file="plots/fitnessModelling_allRegimes_justMAE.pdf", 
       plotFitnessModelling_allRegimes_justMAE, 
       width=8, height=2.75, device="pdf")


