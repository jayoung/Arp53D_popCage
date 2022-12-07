# Arp53D_popCage
population cage modelling for Schroeder et al, 2021
See manuscript for full details.

## Citation:

An actin-related protein that is most highly expressed in Drosophila testes is critical for embryonic development

Courtney M Schroeder, Sarah A Tomlin, Isabel Mejia Natividad, John R Valenzuela, Janet M Young, Harmit S Malik

Division of Basic Sciences & Howard Hughes Medical Institute, Fred Hutchinson Cancer Research Center, Seattle,
United States

eLife 2021;10:e71279 https://elifesciences.org/articles/71279


## Files:

`Arp53D_populationCage_functions.R`  contains all the functions that we'll use for modelling

`Arp53D_populationCage_doModelling.R`   sets up a variety of selection regimes and coefficients, and runs the modelling on each for 35 generations. Takes a few minutes to run. Compares each model with real data using MAE (mean absolute error), chooses the models with lowest MAE (three best models, one for each way of determining heterozygote fitness). Saves results for later plotting.

`Arp53D_populationCage_plotModellingResults.R` makes the plots we're showing in the paper.

## Other resources:

This nice [web tool](https://www.radford.edu/~rsheehy/Gen_flash/popgen/) by Bob Sheehy (Radford University, Virginia) can also perform fitness modelling. We used it to cross-check some of our results.

## Methods (from the paper):

### Population cage setup
The isogenized Arp53D-KO line in the w1118 background was used due to ease of DsRed detection
in the eye (as opposed to the ocelli in the Oregon-R Arp53D-KO background). Virgin females and
males were collected from the w1118 fly line and the Arp53D-KO fly line isogenized in the w1118
background. Crosses with 50 Arp53D-KO females, 25 Arp53D-KO males, and 25 w1118 males were
setup in bottles with three replicates. Crosses were passaged every 2 weeks at room temperature.
At each passage, 50 females and 50 males were randomly collected without fluorescence detection
and without selection based on virgin status. These 100 flies were placed in a fresh bottle, and the
remaining progeny were frozen for subsequent detection of DsRed fluorescence. After 1 week of
laying before the next generation hatched, the 100 flies were removed and frozen to include in the
previous generation’s quantification.

### Modeling selection coefficients
In order to gain insight into the fitness differences between Arp53D-KO and WT flies, we mod-
eled the population cage experiment (Figure 6—figure supplement 1) in silico, simulating experi-
mental evolution using a large number of different fitness parameters (https://github.com/jayoung/
Arp53D_popCage; Young, 2021; copy archived at swh:1:rev:52ff682daab06ba677f43a49de6f5b-
d8a0c54a62). Our modeling assumes a freely mating population of infinite size. We defined fitness
coefficients for each genotype (FWT, Fhet, FKO) relative to WT homozygous flies (FWT = 1). We
explored fitness coefficients for KO homozygous flies (FKO) that ranged between 0.4 and 1 in incre-
ments of 0.001. We explored three possibilities for heterozygote fitness, where fitness matched
either WT (Fhet = FWT), or KO homozygotes (Fhet = FKO), or was exactly intermediate in fitness
between WT and KO homozygotes (Fhet = (FWT + FKO)/2). We seeded all models using the same
genotype combinations as the actual experiment (100% KO homozygous females, and a 50:50 mix
of WT homozygous and KO homozygous males). At each generation, we calculated the fraction of
randomly selected mating pairs that represented each possible genotype combination (Pmat x pat).
For each combination of mating pair genotypes, we used Mendelian segregation to determine the
fraction of offspring genotypes (OWT, Ohet, OKO). To obtain the overall fraction of progeny geno-
types from all parental genotype combinations, we summed the product of those frequencies (P O)
for all mating pair combinations. After obtaining initial progeny genotype frequencies in each gener-
ation, we applied fitness coefficients, multiplying the genotype frequencies by FWT, Fhet, FKO, and re-
normalizing genotype frequencies to sum to 1. This strategy oversimplifies the true biology as it
applies fitness coefficients only to individual genotypes at each generation, regardless of parental
genotypes that we know have strong effects. We iterated these steps over 35 generations and
recorded genotype frequencies at each generation. In order to determine which model best fit the
data, we calculated the mean absolute error (MAE) for each model (by subtracting the modeled
value at the corresponding generation from each real datapoint, taking the absolute value, and then
calculating the mean) and selected the model that minimized MAE.