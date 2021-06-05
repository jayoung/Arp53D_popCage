# Arp53D_popCage
population cage modelling for Schroeder et al, 2021

See manuscript for full details.

Files:

`Arp53D_populationCage_functions.R`  contains all the functions that we'll use for modelling


`Arp53D_populationCage_doModelling.R`   sets up a variety of selection regimes and coefficients, and runs the modelling on each for 35 generations. Takes a few minutes to run. Compares each model with real data using MAE (mean absolute error), chooses the models with lowest MAE (three best models, one for each way of determining heterozygote fitness). Saves results for later plotting.


`Arp53D_populationCage_plotModellingResults.R` makes the plots we're showing in the paper.

