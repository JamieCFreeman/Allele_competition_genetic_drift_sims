# Scott_et-al-2020_Genetic_drift_sims

Based on script written by Melissa Hardstone for:
Hardstone, M. C., B. P. Lazzaro, and J. G.  Scott. 2009. 'The effect of three environmental conditions on the fitness 
  of cytochrome P450 monooxygenase-mediated permethrin resistance in Culex pipiens quinquefasciatus', BMC Evol. Biol., 9: 42.

Intended to test likelihood that allele frequency changes in a population over time could be explained by genetic drift. 

Simulates random mating in a population with four alleles based on supplied starting allele frequencies. 

The simulation assumes a panmictic population with N diploid individuals with allele frequencies defined by the 
observed allele frequencies at the start and end of the interval. For each generation, alleles from the population are sampled with 
replacement to determine the frequency distribution for the next generation. This process is repeated for the number of generations 
defined. p-valuesare obtained by dividing the number of simulations that resulted in an allele frequency 
change as great as or greater than the observed allele frequency at the end of the interval of interest by the total number of 
simulations run. 

Data from the paper is provided in Four_allele_comp_freqs.txt as a working example.

