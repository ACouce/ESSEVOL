  x11(width = 7/2, height = 4, pointsize = 14, canvas = "white")
  plot(c(1:t),x_t[1:t,1],type='o',col='blue',ylim=c(1e-7,1),xlim=c(0,t),log='y',xlab="Time (days)", ylab="Frequency", cex=0.75) #wt (understood as motile-specialist)
  points(c(1:t),yl_t[1:t,1],type='o',col='red', cex=0.75) #non-motile specialist
  points(c(1:t),ye_t[1:t,1],type='o',col='green', cex=0.75) #ess 
  points(c(1:t),yi_t[1:t,1],type='o',col='cyan', cex=0.75) #int generalist
  points(c(1:t),yb_t[1:t,1],type='o',col='magenta', cex=0.75) # bet-hedger1
  points(c(1:t),yb_t[1:t,2],type='o',col='grey', cex=0.75)  # bet-hedger2
