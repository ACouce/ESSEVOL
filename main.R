# .. This program simulates the serial passage of an asexual population in a fluctuating environment in which specialist and generalist mutants can emerge (see README).
# .. This project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No 750129 https://cordis.europa.eu/project/rcn/209320/factsheet/es

  #rm(list=ls())

# .. start the stopwatch
  clock <- Sys.time()
  
  replicas=10
  data=matrix(0,1,replicas)

# .. PARAMETERS
  # .. population sizes
  Nmax0=9.2e+09 # no much diff with 2.3e+10 (14 vs 13 days, no extinctions in neither)
  Nmax1=6.4e+09
  Btlnck0=0.005*Nmax0 #1:200
  Btlnck1=0.000179*Nmax1 #1:5500
  phendel=1
  e<- 2.7182

  # .. mutation rates 
  mu_loss=1e-6*1
  mu_ess=2.2e-10*1 #2.2 × 10−10 Lee 2012 PNAS
  mu_inter=1e-6*1
  mu_bet=1e-9*1
  mu_switch=1e-1

  # .. fitness effects

  LR=1.285 #empirical estimates
  LP=1.41 #minimum 1.41 #average 2.096
  SR=0.43 #average 0.336
  SP=0.557

  W=LP # s=w/ln(2)-1 #link to pop genet
  s_ben=e^(W*log(2))-2 
  W=SR
  s_cost=e^(W*log(2))-2     
  i_eff=0.75

  # .. environmental pattern
  tmax<-100
  medium<-rep(c(0,0,0,1),tmax)

  # .. genotypic matrices
  n_ben=1
  n_del=1
  rx=matrix(nrow = n_ben, ncol = n_del)
  ryl=matrix(nrow = n_ben, ncol = n_del)
  rye=matrix(nrow = n_ben, ncol = n_del)
  ryi=matrix(nrow = n_ben, ncol = n_del)
  ryb=matrix(nrow = n_ben, ncol = 2)

  # .. VARIABLES
  x=matrix(nrow = n_ben, ncol = n_del)
  yl=matrix(nrow = n_ben, ncol = n_del)
  ye=matrix(nrow = n_ben, ncol = n_del)
  yi=matrix(nrow = n_ben, ncol = n_del)
  yb=matrix(nrow = n_ben, ncol = 2)

  ratiox=matrix(nrow = n_ben, ncol = n_del)
  x_t=matrix(nrow = tmax+1,  ncol = 1)

  ratioyl=matrix(nrow = n_ben, ncol = n_del)
  yl_t=matrix(nrow = tmax+1,  ncol = 1) 

  ratioye=matrix(nrow = n_ben, ncol = n_del)
  ye_t=matrix(nrow = tmax+1,  ncol = 1)

  ratioyi=matrix(nrow = n_ben, ncol = n_del)
  yi_t=matrix(nrow = tmax+1,  ncol = 1)

  ratioyb=matrix(nrow = n_ben, ncol = 2)
  yb_t=matrix(nrow = tmax+1,  ncol = 2)

  for (replica in 1:replicas) {
  # .. ALGORITHM
  t=0
  g=0
  x[]=0
  x[1,1]=Btlnck0
  yl[]=0
  ye[]=0
  yi[]=0
  yb[]=0
  yb[1,1]=0
  ratiox=x/(sum(x)+sum(yl)+sum(ye)+sum(yi)+sum(yb))
  ratioyl=yl/(sum(x)+sum(yl)+sum(ye)+sum(yi)+sum(yb))
  ratioye=ye/(sum(x)+sum(yl)+sum(ye)+sum(yi)+sum(yb))
  ratioyi=yi/(sum(x)+sum(yl)+sum(ye)+sum(yi)+sum(yb))
  ratioyb=yb/(sum(x)+sum(yl)+sum(ye)+sum(yi)+sum(yb))

  source("serial_culture.R")
  if (replica%in%seq(1,100,5)) {print(c(replica,t))}
  data[replica]<-t

  # .. GRAPHIC OUTPUT

  source("plotter.R")
  

}
  hist(data, xlim=c(1,30), breaks=15, col="grey")
  print(c(mean(data[data<tmax]), sum(data==replica)))
# .. stop the stopwatch1
  print(Sys.time()-clock)


