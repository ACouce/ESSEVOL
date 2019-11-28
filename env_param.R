
  if (environment==0) { #shaken
  Nmax=Nmax0

	  for (i in 1:n_ben) {
	      for (j in 1:n_del) { 
		   rx[i,j]=2
		   ryl[i,j]=2+s_ben
		   rye[i,j]=2+s_ben
		   ryi[i,j]=2+s_ben*i_eff
		   ryb[]<-c(2,2+s_ben)
	     			 }
			  }
  }  else if (environment==1) { #structured
  Nmax=Nmax1
	for (i in 1:n_ben) {
		for (j in 1:n_del) { 
		   rx[i,j]=2
		   ryl[i,j]=2+s_cost
		   rye[i,j]=2
		   ryi[i,j]=2+s_cost*i_eff
		   ryb[]<-c(2,2+s_cost)
		     		}
			}
  }

#print(c("environment=",environment))
#print(c(rx, ryl, rye, ryi, ryb))
