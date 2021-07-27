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

