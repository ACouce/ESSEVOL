
 


	if (g<phendel) { #generations to express ESS phenotype 

		  if (environment==0) { #now is shaken, but phenotype is motile
			   rye[1,1]=2

		  }  else if (environment==1) { #now is structured, but phenotype is non-motile
			   rye[1,1]=2+s_cost

		  }

	     } else {

		  if (environment==0) { #shaken
			   rye[1,1]=2+s_ben

		  }  else if (environment==1) { #structured
			   rye[1,1]=2

		  }

	}

#print('change!')
#print(c(t,g,medium[t], medium[t-1], environment, rye))
