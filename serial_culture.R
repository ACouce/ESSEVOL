# .. This routine simulates a single serial culture experiment. It is intended to be 
# .. called by the script 'main.R'

#while ((sum(ratiox[2,]+ratioy[2,]) < 1) && (sum(ratioys[2,]) < 1)) {
while ((t<tmax) & (sum(ratioye[]) < 0.5)) {

  t=t+1
  environment<-medium[t]
  source("env_param.R")
  g=0	
  # .. this loop represents one day in a flask 
  while (sum(x)+sum(yl)+sum(ye)+sum(yi)+sum(yb) < Nmax) {
    #print(c(g, sum(x)+sum(yl)+sum(ye)+sum(yi)+sum(yb)))
    if (t>1) { 
        if (medium[t] != medium[t-1]) {source("phen_delay.R")} #only if we just changed environment
    }    
    g=g+1
    x[]=rx[]*x
    yl[]=ryl[]*yl
    ye[]=rye[]*ye
    yi[]=ryi[]*yi
    yb[]=ryb[]*yb
    
    for (i in 1:n_ben) {
	for (j in 1:n_del) {

	#.. allocation of mutations
	mloss<- rpois(1,lambda=x[i,j]*mu_loss)
	mess<- rpois(1,lambda=x[i,j]*mu_ess)
	minter<- rpois(1,lambda=x[i,j]*mu_inter)
	mbet<- rpois(1,lambda=x[i,j]*mu_bet)

	#.. bet-hedging switches

    if (yb[i,1]*mu_switch>.Machine$integer.max) {
        mswitchF<-round(rnorm(1,mean=yb[i,1]*mu_switch,sd=sqrt(yb[i,1]*mu_switch)))
    } else {
        mswitchF<- rpois(1,lambda=yb[i,1]*mu_switch)
        }

    if (yb[i,2]*mu_switch>.Machine$integer.max) {
        mswitchR<-round(rnorm(1,mean=yb[i,2]*mu_switch,sd=sqrt(yb[i,2]*mu_switch)))
    } else {
        mswitchR<- rpois(1,lambda=yb[i,2]*mu_switch)
        }

   # print('potato')
  #  print(c(yb[i,1], yb[i,2], mswitchF, mswitchR, ryb))
	#.. this if is for the boundaries of the genotypic matrix
	a=1
	b=1

	if (i==n_ben) { a=0 }
	if (j==n_del) { b=0 }

	x[i,j]=x[i,j]-mloss-mess-minter-mbet


	yl[i,j]=yl[i,j]+mloss
	ye[i,j]=ye[i,j]+mess 
	yi[i,j]=yi[i,j]+minter

	yb[i,1]=yb[i,1]+mbet+mswitchR-mswitchF  
    yb[i,2]=yb[i,2]+mswitchF-mswitchR      
 
   #print(c(t,g,yb, yb[i,1]+mbet+mswitchR-mswitchF  , yb[i,2]+mswitchF-mswitchR   ))
	
	x[x<0]=0
	yl[yl<0]=0
	ye[ye<0]=0
	yi[yi<0]=0
	yb[yb<0]=0
   # print(yb)
	}
    }
  }
    ratiox=x/(sum(x)+sum(yl)+sum(ye)+sum(yi)+sum(yb))
    ratioyl=yl/(sum(x)+sum(yl)+sum(ye)+sum(yi)+sum(yb))
    ratioye=ye/(sum(x)+sum(yl)+sum(ye)+sum(yi)+sum(yb))
    ratioyi=yi/(sum(x)+sum(yl)+sum(ye)+sum(yi)+sum(yb))
    ratioyb1=yb[1,1]/(sum(x)+sum(yl)+sum(ye)+sum(yi)+sum(yb))
    ratioyb2=yb[1,2]/(sum(x)+sum(yl)+sum(ye)+sum(yi)+sum(yb))

	x_t[t,1]<-(sum(ratiox[1,]))
	yl_t[t,1]<-(sum(ratioyl[1,]))
	ye_t[t,1]<-(sum(ratioye[1,]))
	yi_t[t,1]<-(sum(ratioyi[1,]))
	yb_t[t,1]<-ratioyb1
	yb_t[t,2]<-ratioyb2

   # .. the bottleneck
    Btlnck=Btlnck0*(1-medium[t+1])+medium[t+1]*(Btlnck1+medium[t+1]-1)
    #print(c(medium[t], Btlnck))
    for (i in 1:n_ben) {
	for (j in 1:n_del) {
     # print(c(sum(x)+sum(yl)+sum(ye)+sum(yi)+sum(yb),Nmax, Btlnck))
	  if ((ratiox[i,j]>0.1)) {
	    x[i,j]=rbinom (1, Btlnck, ratiox[i,j])
	  }
	  else if ((ratiox[i,j]>0)) {

	    x[i,j]=rpois (1, lambda=Btlnck*ratiox[i,j])
	  } 
	  else {
            x[i,j]=0
	  }

	  if ((ratioyl[i,j]>0.1)) {
	    yl[i,j]=rbinom (1, Btlnck, ratioyl[i,j])
	  }
	  else if ((ratioyl[i,j]>0)) {

	    yl[i,j]=rpois (1, lambda=Btlnck*ratioyl[i,j])
	  } 
	  else {
            yl[i,j]=0
	  }

	  if ((ratioye[i,j]>0.1)) {
	    ye[i,j]=rbinom (1, Btlnck, ratioye[i,j])
	  }
	  else if ((ratioye[i,j]>0)) {

	    ye[i,j]=rpois (1, lambda=Btlnck*ratioye[i,j])
	  } 
	  else {
            ye[i,j]=0
	  }

	  if ((ratioyi[i,j]>0.1)) {
	    yi[i,j]=rbinom (1, Btlnck, ratioyi[i,j])
	  }
	  else if ((ratioyi[i,j]>0)) {

	    yi[i,j]=rpois (1, lambda=Btlnck*ratioyi[i,j])
	  } 
	  else {
            yi[i,j]=0
	  }
	  if ((ratioyb1>0.1)) {
	    yb[i,1]=rbinom (1, Btlnck, ratioyb1)
	  }
	  else if ((ratioyb1>0)) {

	    yb[i,1]=rpois (1, lambda=Btlnck*ratioyb1)
	  } 
	  else {
            yb[i,1]=0
	  }
	  if ((ratioyb2>0.1)) {
	    yb[i,2]=rbinom (1, Btlnck, ratioyb2)
	  }
	  else if ((ratioyb2>0)) {

	    yb[i,2]=rpois (1, lambda=Btlnck*ratioyb2)
	  } 
	  else {
            yb[i,2]=0
	  }
	}
     }

}

