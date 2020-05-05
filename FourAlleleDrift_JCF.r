
# This function is intended to simulate random drift of 4 alleles in a population, given starting and ending
# frequencies for each allele, by resampling population of given allele frequency.
# Revision of a script by Melissa Hardstone (FourAlleleDrift.r), reviewed and minor edits JCF 2018.
# The MH copy of FourAlleleDrift.r script HAS AN ERROR- in the reporting of results for A2 and can give incorrect results.
# What I changed:
    #
    # Results were previously only kept for current simulation, I routed them to a larger matrix that stores 
    #      results for all sims during function. If raw simulation data is desired, default "raw_data=FALSE" can 
    #      be changed to true and all data will be returned. (Good for making sure everything is working as intended, 
    #      and plotting if desired.)
    # Also quieted down the function a bit, there were unnecessary print commands everywhere.
    # Added comments for clarity
    

# last edit: 2019-1-8 JCF

# R sessionInfo of last edit.
# > sessionInfo()
# R version 3.5.1 (2018-07-02)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows >= 8 x64 (build 9200)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# loaded via a namespace (and not attached):
#   [1] compiler_3.5.1 tools_3.5.1   



#
# Input:
#
# SIMS is # simulations desired
# gens is # generations between starting & ending allele frequencies
# size_of_pop is # of diploid individuals in the population
# A1-A4 are starting allele frequencies   
# A1-A4 are starting allele frequencies

#
# Output:
#
# Vector of 4 p-values (1 for each allele), calculated as the number of simulations where the allele frequency change is as extreme as the
#         observed allele frequency change. Null hypothesis here is that allele frequency change could be due to drift.
# raw_data is set to FALSE as default, but if set to TRUE function will return a matrix 


FourAlleleDrift<-function(SIMS, gens, size_of_pop, A1, A2, A3, A4, end_A1, end_A2, end_A3, end_A4, raw_data=FALSE)
{
  
  # Set up 4 matrices to track allele frequency over gens (dimensions #gens X #sims)
  # Fill first row with the starting allele frequency
	resultA1<-matrix(NA, nrow=gens, ncol=SIMS)                    # Declare results matrix
	resultA1[1,]=A1                                               # Fill first row with the starting allele frequency
	resultA2<-matrix(NA, nrow=gens, ncol=SIMS)
	resultA2[1,]=A2
	resultA3<-matrix(NA, nrow=gens, ncol=SIMS)
	resultA3[1,]=A3
  resultA4<-matrix(NA, nrow=gens, ncol=SIMS)
  resultA4[1,]=A4

  
  # Now will loop through SIMS (j) & generations (i)
  # 
	for (j in 1:SIMS)                                              # For each of j simulations
	{

	  
    gen<-gens                                                    # no. of generations
		pop_size<-size_of_pop                                        # size_of_pop
		i=1                                                          # Initialize for generation

		# starting with the inputed allele frequencies
		A1_allele<-A1                                                # A1 start freq
		A2_allele<-A2                                                # A2 start freq
		A3_allele<-A3                                                # A3 start freq
		A4_allele<-A4                                                # A4 start freq

		# Create a population with the starting allele frequencies (*2 for diploid)
		pop_alleles<-c( rep('A1',2*pop_size*A1_allele), rep('A2',2*pop_size*A2_allele),
		                rep('A3',2*pop_size*A3_allele), rep('A4',2*pop_size*A4_allele) )

		# take a sample with replacement from the whole population to make next generation (panmictic, random mating)
		pop_alleles<-sample( pop_alleles, 2*pop_size, replace=TRUE )

		# Add the simulated allele frequency to results matrices: find # of entries in pop_alleles for each allele & divide by 2*pop_size
		resultA1[i+1,j]<-( length(pop_alleles[pop_alleles=='A1']) )/(2*pop_size)   # freq of A1 allele, add to results matrix
		resultA2[i+1,j]<-( length(pop_alleles[pop_alleles=='A2']) )/(2*pop_size)   # freq of A2 allele
		resultA3[i+1,j]<-( length(pop_alleles[pop_alleles=='A3']) )/(2*pop_size)   # freq of A3 allele
		resultA4[i+1,j]<-( length(pop_alleles[pop_alleles=='A4']) )/(2*pop_size)   # freq of A4 allele
	
		
		
		
    # each loop, allele freq ending from previous loop used as the start allele freq
		for(i in 2:gen-1)
		{
		  A1_allele<-resultA1[(i+1)-1,j]                                           # uses A1 freq from last generation
			A2_allele<-resultA2[(i+1)-1,j]                                           # uses new A2 freq
			A3_allele<-resultA3[(i+1)-1,j]                                           # uses new A3 freq
			A4_allele<-resultA4[(i+1)-1,j]                                           # uses new A4 freq

      # Simulate population with new allele frequencies
      pop_alleles<-c(rep('A1',2*pop_size*A1_allele),rep('A2',2*pop_size*A2_allele),rep('A3',2*pop_size*A3_allele),rep('A4',2*pop_size*A4_allele));

			# take a sample from the whole population to make next generation
      pop_alleles<-sample(pop_alleles,2*pop_size,replace=TRUE)
      resultA1[i+1,j]<-(length(pop_alleles[pop_alleles=='A1']))/(2*pop_size)   # freq of A1 allele, add to results matrix
      resultA2[i+1,j]<-(length(pop_alleles[pop_alleles=='A2']))/(2*pop_size)   # freq of A2 allele
      resultA3[i+1,j]<-(length(pop_alleles[pop_alleles=='A3']))/(2*pop_size)   # freq of A3 allele
      resultA4[i+1,j]<-(length(pop_alleles[pop_alleles=='A4']))/(2*pop_size)   # freq of A4 allele

	  }

	} 
  
  # Now have full matrix of results for each allele.
  # Now calculate p-values from simulation results.

  # count times A1 allele freq at generation of interest is equal to or more extreme than obs freq
    # need to designate which direction of the extreme values is wanted
    if(end_A1 < A1)
    {
		  A1_p<-length( which( resultA1[gen,] <= end_A1 ) )/SIMS;
		  #cat("probability of getting end A1 freq or smaller is", A1_p, '\n');
	 }else{
	   A1_p<-length( which( resultA1[gen,] >= end_A1 ) )/SIMS;
		  #cat("probability of getting end A1 freq or larger is", A1_p, '\n');
	 }

 	# Count times A2 allele freq at generation of interest is equal to or more extreme than obs freq
	# Need to designate which direction of the extreme values is wanted
  if(end_A2 < A2)
  {
   A2_p<-length( which( resultA2[gen,] <= end_A2 ) )/SIMS;
    # cat("probability of getting end A2 freq or smaller is", A2_p, '\n');
  }else{
    A2_p<-length( which( resultA2[gen,] >= end_A2 ) )/SIMS;
    # cat("probability of getting end A2 freq or larger is", A2_p, '\n');
  }

	# count times A3 allele freq at generation of interest is equal to or more extreme than obs freq
	# need to designate which direction of the extreme values is wanted
  if(end_A3 < A3)
  {
    A3_p<-length( which( resultA3[gen,] <= end_A3 ) )/SIMS;
    #cat("probability of getting end A3 freq or smaller is", A3_p, '\n');
  }else{
    A3_p<-length( which( resultA3[gen,] >= end_A3 ) )/SIMS;
    #cat("probability of getting end A3 freq or larger is", A3_p, '\n');
  }
	
	# count times A4 allele freq at generation of interest is equal to or more extreme than obs freq
	# need to designate which direction of the extreme values is wanted
  if(end_A4 < A4)
  {
    A4_p<-length( which( resultA4[gen,] <= end_A4 ) )/SIMS;
    #cat("probability of getting end A4 freq or smaller is", A4_p, '\n');
  }else{
    A4_p<-length( which( resultA4[gen,] >= end_A4 ) )/SIMS;
    #cat("probability of getting end A4 freq or larger is", A4_p, '\n');
  }
  
  
  # Function by default returns a vector of 4 p-values, 1 for each allele. Null hypothesis here is that allele frequency change could be due to drift.
  # If raw_data==TRUE, will also return a matrix of simulation results (intended for plotting). 
  
  # Create vector of p-values
  p.vals=c(A1_p,A2_p,A3_p,A4_p)
  
  
  if( raw_data==TRUE ){
    print( list("A1_results"=resultA1, "A2_results"=resultA2, "A3_results"=resultA3, "A4_results"=resultA4, p.vals) )
  } else { print(p.vals) }
  
  
  
}





