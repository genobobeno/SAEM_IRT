PlotErrorEstimates<-function(ErrorStats) {

  stats1D<-ErrorStats
  nf <- layout(matrix(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,8,8,8,8,8,9,9,9,9,9,9,9),14,2), c(18,24), 
               c(4,3,3,3,3,3,3,3,3,3,3,3,5,4,3,3,3,3,3,3,3,3,3,3,3,5), TRUE)
  #layout.show(nf)

  YR<-range(c(#stats1D$A1_ERR_CLT_Ann,
    stats1D$A1_ERR_CLT_Burn,
    #stats1D$A1_ERR_CLT_Emp,
    stats1D$A1_ERR_LMI_Iter,
    stats1D$A1_ERR_LMI_Iter_Emp_2PL,
    stats1D$A1_ERR_LMI_Iter_Emp_2PNO,
    stats1D$A1_ERR_Simple_2PL,
    stats1D$A1_ERR_Simple_2PNO,
    #stats1D$A1_ERR_MCMC_Ann,
    stats1D$A1_ERR_MCMC_Burn #,
    #stats1D$A1_ERR_MCMC_Emp
  )/stats1D$A1_RMSE)
  
  txt.labels<-c(#expression("CLT"["Ann"]),
    expression("CLT"["B"]),
    #expression("CLT"["Emp"]),
    #expression("MCMC"["Ann"]),
    expression("MCMC"["B"]),
    #expression("MCMC"["Emp"]),
    expression("ICE"),
    expression("IPCE"["L"]),
    expression("IPCE"["O"]),
    expression("SPCE"["L"]),
    expression("SPCE"["O"]))
  
  ticks <- pretty(YR)[-c(length(pretty(YR))-1:0)]
  labels <- format(ticks, big.mark=",", scientific=FALSE)
  m.b = 0.25
  par(mar=c(m.b,4,2,m.b),cex=0.9)
  # for (p in 1:ncol(gen.xi)) {
  p=1
  ifelse(p==ncol(gen.xi),pt<-"B",pt<-paste("A",p,sep=""))
  # plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
  #      main=expression("Slope Error Estimates : 2PNO"),ylim=YR,
  #      xlab=NA,ylab=NA,xaxt="n",yaxt="n")
  # axis(2,at=ticks,labels=labels,las=2)
  # abline(h=1,lwd=2,lty=2)
  # points(gen.xi[,p],stats1D$A1_ERR_CLT_Ann/stats1D$A1_RMSE,pch=0,col=1)
  # text(0.6,0.9*YR[2],txt.labels[1])
  
  #par(mar=c(m.b,4,m.b,m.b))
  plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
       main=expression("Slope Error Estimates : 2PNO"),
       ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
  axis(2,at=ticks,labels=labels,las=2)
  abline(h=1,lwd=2,lty=2)
  points(gen.xi[,p],stats1D$A1_ERR_CLT_Burn/stats1D$A1_RMSE,pch=1,col=1)
  text(0.6,0.9*YR[2],txt.labels[1])
  
  par(mar=c(m.b,4,m.b,m.b))
  # plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
  #      main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
  # axis(2,at=ticks,labels=labels,las=2)
  # abline(h=1,lwd=2,lty=2)
  # points(gen.xi[,p],stats1D$A1_ERR_CLT_Emp/stats1D$A1_RMSE,pch=2,col=1)
  # text(0.6,0.9*YR[2],txt.labels[3])
  
  # plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
  #      main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
  # axis(2,at=ticks,labels=labels,las=2)
  # abline(h=1,lwd=2,lty=2)
  # points(gen.xi[,p],stats1D$A1_ERR_MCMC_Ann/stats1D$A1_RMSE,pch=8,col=1)
  # text(0.6,0.9*YR[2],txt.labels[4])
  
  plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
       main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
  axis(2,at=ticks,labels=labels,las=2)
  abline(h=1,lwd=2,lty=2)
  points(gen.xi[,p],stats1D$A1_ERR_MCMC_Burn/stats1D$A1_RMSE,pch=9,col=1)
  text(0.6,0.9*YR[2],txt.labels[2])
  
  # plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
  #      main=NA,ylim=YR,xlab=NA,ylab=NA,yaxt="n",xaxt="n")
  # axis(2,at=ticks,labels=labels,las=2)
  # abline(h=1,lwd=2,lty=2)
  # points(gen.xi[,p],stats1D$A1_ERR_MCMC_Emp/stats1D$A1_RMSE,pch=10,col=1)
  # text(0.6,0.9*YR[2],txt.labels[6])
  
  plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
       main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
  axis(2,at=ticks,labels=labels,las=2)
  abline(h=1,lwd=2,lty=2)
  points(gen.xi[,p],stats1D$A1_ERR_LMI_Iter/stats1D$A1_RMSE,pch=3,col=1)
  text(0.6,0.9*YR[2],txt.labels[3])
  
  plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
       main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
  axis(2,at=ticks,labels=labels,las=2)
  abline(h=1,lwd=2,lty=2)
  points(gen.xi[,p],stats1D$A1_ERR_LMI_Iter_Emp_2PL/stats1D$A1_RMSE,pch=4,col=1)
  text(0.6,0.9*YR[2],txt.labels[4])
  mtext(text = "Ratio of Error Approximation to RMSE",side = 2,line = 2.5,cex=0.7,adj=0)
  
  plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
       main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
  axis(2,at=ticks,labels=labels,las=2)
  abline(h=1,lwd=2,lty=2)
  points(gen.xi[,p],stats1D$A1_ERR_LMI_Iter_Emp_2PNO/stats1D$A1_RMSE,pch=5,col=1)
  text(0.6,0.9*YR[2],txt.labels[5])
  
  plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
       main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
  axis(2,at=ticks,labels=labels,las=2)
  abline(h=1,lwd=2,lty=2)
  points(gen.xi[,p],stats1D$A1_ERR_Simple_2PL/stats1D$A1_RMSE,pch=6,col=1)
  text(0.6,0.9*YR[2],txt.labels[6])
  
  par(mar=c(4+0.5,4,m.b,m.b))
  plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
       main=NA,ylim=YR,xlab="Generated Slope (1D, 2PNO)",ylab=NA,yaxt="n")
  axis(2,at=ticks,labels=labels,las=2)
  abline(h=1,lwd=2,lty=2)
  points(gen.xi[,p],stats1D$A1_ERR_Simple_2PNO/stats1D$A1_RMSE,pch=7,col=1)
  text(0.6,0.9*YR[2],txt.labels[7])
  
  
  # fitA1.CLT_Ann<-loess(A1_ERR_CLT_Ann/A1_RMSE~A1,data=stats1D,span=0.5)
  # fitA1.CLT_Ann
  # sd(fitA1.CLT_Ann$residuals)
  #fitA1o<-lm(A1_sd_LMIo/A1_sd~A1,data=stats1D)
  fitA1.CLT_Burn<-loess(A1_ERR_CLT_Burn/A1_RMSE~A1,data=stats1D,span=0.5)
  fitA1.CLT_Burn
  sd(fitA1.CLT_Burn$residuals)
  #fitA1MC<-lm(A1_sd_MCE/A1_sd~A1,data=stats1D)
  # fitA1.CLT_Emp<-loess(A1_ERR_CLT_Emp/A1_RMSE~A1,data=stats1D,span=0.5)
  # fitA1.CLT_Emp
  # sd(fitA1.CLT_Emp$residuals)
  #fitA1LM<-lm(A1_sd_LMI/A1_sd~A1,data=stats1D)
  fitA1.LMI_Iter<-loess(A1_ERR_LMI_Iter/A1_RMSE~A1,data=stats1D,span=0.5)
  fitA1.LMI_Iter
  sd(fitA1.LMI_Iter$residuals)
  fitA1.LMI_Iter_Emp_2PL<-loess(A1_ERR_LMI_Iter_Emp_2PL/A1_RMSE~A1,data=stats1D,span=0.5)
  fitA1.LMI_Iter_Emp_2PL
  sd(fitA1.LMI_Iter_Emp_2PL$residuals)
  #fitA1o<-lm(A1_sd_LMIo/A1_sd~A1,data=stats1D)
  fitA1.LMI_Iter_Emp_2PNO<-loess(A1_ERR_LMI_Iter_Emp_2PNO/A1_RMSE~A1,data=stats1D,span=0.5)
  fitA1.LMI_Iter_Emp_2PNO
  sd(fitA1.LMI_Iter_Emp_2PNO$residuals)
  #fitA1MC<-lm(A1_sd_MCE/A1_sd~A1,data=stats1D)
  fitA1.Simple_2PL<-loess(A1_ERR_Simple_2PL/A1_RMSE~A1,data=stats1D,span=0.5)
  fitA1.Simple_2PL
  sd(fitA1.Simple_2PL$residuals)
  #fitA1LM<-lm(A1_sd_LMI/A1_sd~A1,data=stats1D)
  fitA1.Simple_2PNO<-loess(A1_ERR_Simple_2PNO/A1_RMSE~A1,data=stats1D,span=0.5)
  fitA1.Simple_2PNO
  sd(fitA1.Simple_2PNO$residuals)
  # fitA1.MCMC_Ann<-loess(A1_ERR_MCMC_Ann/A1_RMSE~A1,data=stats1D,span=0.5)
  # fitA1.MCMC_Ann
  # sd(fitA1.MCMC_Ann$residuals)
  #fitA1o<-lm(A1_sd_LMIo/A1_sd~A1,data=stats1D)
  fitA1.MCMC_Burn<-loess(A1_ERR_MCMC_Burn/A1_RMSE~A1,data=stats1D,span=0.5)
  fitA1.MCMC_Burn
  sd(fitA1.MCMC_Burn$residuals)
  #fitA1MC<-lm(A1_sd_MCE/A1_sd~A1,data=stats1D)
  # fitA1.MCMC_Emp<-loess(A1_ERR_MCMC_Emp/A1_RMSE~A1,data=stats1D,span=0.5)
  # fitA1.MCMC_Emp
  # sd(fitA1.MCMC_Emp$residuals)
  ######plot2
  par(mar=c(4,6,4,m.b))
  seqA1<-seq(range(stats1D$A1)[1],range(stats1D$A1)[2],length.out = 50)
  plot(stats1D$A1,stats1D$A1_ERR_MCMC_Burn/stats1D$A1_RMSE,main=expression("Loess fit to"~hat(sigma)/sigma["RMSE"]),
       ylab="Ratio of Error Approximation to RMSE",xlab="Generated Slope (1D, 2PNO)",ylim=YR,type="n")
  #lines(seqA1,predict(fitA1.CLT_Ann,newdata = data.frame(A1=seqA1)),lty=1)
  lines(seqA1,predict(fitA1.CLT_Burn,newdata = data.frame(A1=seqA1)),lty=2,col=2)
  #lines(seqA1,predict(fitA1.CLT_Emp,newdata = data.frame(A1=seqA1)),lty=3,col=3)
  #lines(seqA1,predict(fitA1.MCMC_Ann,newdata = data.frame(A1=seqA1)),lty=4,col=4,lwd=2)
  lines(seqA1,predict(fitA1.MCMC_Burn,newdata = data.frame(A1=seqA1)),lty=5,col=6,lwd=2)
  #lines(seqA1,predict(fitA1.MCMC_Emp,newdata = data.frame(A1=seqA1)),lty=6,col=1)
  legend("topright",c(#expression("CLT"["Ann"]),
    expression("CLT"["B"]),
    #expression("CLT"["Emp"]),
    #expression("MCMC"["Ann"]),
    expression("MCMC"["B"])#,
    #expression("MCMC"["Emp"])
  ),#inset = 0.05,lty=1:6,col=c(1:4,6,1),lwd=c(1,1,1,2,2,1))
  inset = 0.05,lty=c(2,5),col=c(2,6),lwd=c(1,2))
  abline(h=1)
  par(mar=c(4+0.5,6,7,m.b))
  plot(stats1D$A1,stats1D$A1_ERR_MCMC_Burn/stats1D$A1_RMSE,main=expression("Loess fit to Hessian"~hat(sigma)/sigma["RMSE"]),
       ylab="Ratio of Error Approximation to RMSE",xlab="Generated Slope (1D, 2PNO)",ylim=YR,type="n")
  lines(seqA1,predict(fitA1.LMI_Iter,newdata = data.frame(A1=seqA1)),lty=1,col=1)
  lines(seqA1,predict(fitA1.LMI_Iter_Emp_2PL,newdata = data.frame(A1=seqA1)),lty=2,col=2)
  lines(seqA1,predict(fitA1.LMI_Iter_Emp_2PNO,newdata = data.frame(A1=seqA1)),lty=3,col=3,lwd=2)
  lines(seqA1,predict(fitA1.Simple_2PL,newdata = data.frame(A1=seqA1)),lty=4,col=4)
  lines(seqA1,predict(fitA1.Simple_2PNO,newdata = data.frame(A1=seqA1)),lty=5,col=6)
  legend("topright",c(expression("ICE"),
                      expression("IPCE"["L"]),
                      expression("IPCE"["O"]),
                      expression("SPCE"["L"]),
                      expression("SPCE"["O"])),inset = 0.05,lty=1:5,col=c(1:4,6),lwd=c(1,1,2,1,1))
  abline(h=1)
  
  #########  Error Plots for B
  p=2
  YR<-range(c(#stats1D$B_ERR_CLT_Ann,
    stats1D$B_ERR_CLT_Burn,
    #stats1D$B_ERR_CLT_Emp,
    stats1D$B_ERR_LMI_Iter,
    stats1D$B_ERR_LMI_Iter_Emp_2PL,
    stats1D$B_ERR_LMI_Iter_Emp_2PNO,
    stats1D$B_ERR_Simple_2PL,
    stats1D$B_ERR_Simple_2PNO,
    #stats1D$B_ERR_MCMC_Ann,
    stats1D$B_ERR_MCMC_Burn#,
    #stats1D$B_ERR_MCMC_Emp
  )/stats1D$B_RMSE)
  
  ticks <- pretty(YR)[-c(length(pretty(YR))-1:0)]
  labels <- format(ticks, big.mark=",", scientific=FALSE)
  m.b = 0.25
  par(mar=c(m.b,4,2,m.b))
  # for (p in 1:ncol(gen.xi)) {
  p=2
  ifelse(p==ncol(gen.xi),pt<-"B",pt<-paste("A",p,sep=""))
  # plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
  #      main=expression("Intercept Error Estimates : 2PNO"),ylim=YR,
  #      xlab=NA,ylab=NA,xaxt="n",yaxt="n")
  # axis(2,at=ticks,labels=labels,las=2)
  # abline(h=1,lwd=2,lty=2)
  # points(gen.xi[,p],stats1D$B_ERR_CLT_Ann/stats1D$B_RMSE,pch=0,col=1)
  # text(-1,0.9*YR[2],txt.labels[1])
  
  # par(mar=c(m.b,4,m.b,m.b))
  plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
       main=expression("Intercept Error Estimates : 2PNO"),
       ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
  axis(2,at=ticks,labels=labels,las=2)
  abline(h=1,lwd=2,lty=2)
  points(gen.xi[,p],stats1D$B_ERR_CLT_Burn/stats1D$B_RMSE,pch=1,col=1)
  text(-1,0.9*YR[2],txt.labels[1])
  
  par(mar=c(m.b,4,m.b,m.b))
  # plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
  #      main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
  # axis(2,at=ticks,labels=labels,las=2)
  # abline(h=1,lwd=2,lty=2)
  # points(gen.xi[,p],stats1D$B_ERR_CLT_Emp/stats1D$B_RMSE,pch=2,col=1)
  # text(-1,0.9*YR[2],txt.labels[3])
  
  # plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
  #      main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
  # axis(2,at=ticks,labels=labels,las=2)
  # abline(h=1,lwd=2,lty=2)
  # points(gen.xi[,p],stats1D$B_ERR_MCMC_Ann/stats1D$B_RMSE,pch=8,col=1)
  # text(-1,0.9*YR[2],txt.labels[4])
  
  plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
       main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
  axis(2,at=ticks,labels=labels,las=2)
  abline(h=1,lwd=2,lty=2)
  points(gen.xi[,p],stats1D$B_ERR_MCMC_Burn/stats1D$B_RMSE,pch=9,col=1)
  text(-1,0.9*YR[2],txt.labels[2])
  
  # plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
  #      main=NA,ylim=YR,xlab=NA,ylab=NA,yaxt="n",xaxt="n")
  # axis(2,at=ticks,labels=labels,las=2)
  # abline(h=1,lwd=2,lty=2)
  # points(gen.xi[,p],stats1D$B_ERR_MCMC_Emp/stats1D$B_RMSE,pch=10,col=1)
  # text(-1,0.9*YR[2],txt.labels[6])
  
  plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
       main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
  axis(2,at=ticks,labels=labels,las=2)
  abline(h=1,lwd=2,lty=2)
  points(gen.xi[,p],stats1D$B_ERR_LMI_Iter/stats1D$B_RMSE,pch=3,col=1)
  text(-1,0.9*YR[2],txt.labels[3])
  
  plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
       main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
  axis(2,at=ticks,labels=labels,las=2)
  abline(h=1,lwd=2,lty=2)
  points(gen.xi[,p],stats1D$B_ERR_LMI_Iter_Emp_2PL/stats1D$B_RMSE,pch=4,col=1)
  text(-1,0.9*YR[2],txt.labels[4])
  mtext(text = "Ratio of Error Approximation to RMSE",side = 2,line = 2.5,cex=0.7,adj=0)
  
  plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
       main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
  axis(2,at=ticks,labels=labels,las=2)
  abline(h=1,lwd=2,lty=2)
  points(gen.xi[,p],stats1D$B_ERR_LMI_Iter_Emp_2PNO/stats1D$B_RMSE,pch=5,col=1)
  text(-1,0.9*YR[2],txt.labels[5])
  
  plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
       main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
  axis(2,at=ticks,labels=labels,las=2)
  abline(h=1,lwd=2,lty=2)
  points(gen.xi[,p],stats1D$B_ERR_Simple_2PL/stats1D$B_RMSE,pch=6,col=1)
  text(-1,0.9*YR[2],txt.labels[6])
  
  par(mar=c(4+0.5,4,m.b,m.b))
  plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
       main=NA,ylim=YR,xlab="Generated Intercept (2PNO)",ylab=NA,yaxt="n")
  axis(2,at=ticks,labels=labels,las=2)
  abline(h=1,lwd=2,lty=2)
  points(gen.xi[,p],stats1D$B_ERR_Simple_2PNO/stats1D$B_RMSE,pch=7,col=1)
  text(-1,0.9*YR[2],txt.labels[7])
  
  
  # fitB.CLT_Ann<-loess(B_ERR_CLT_Ann/B_RMSE~B,data=stats1D,span=0.5)
  # fitB.CLT_Ann
  # sd(fitB.CLT_Ann$residuals)
  #fitBo<-lm(B_sd_LMIo/B_sd~B,data=stats1D)
  fitB.CLT_Burn<-loess(B_ERR_CLT_Burn/B_RMSE~B,data=stats1D,span=0.5)
  fitB.CLT_Burn
  sd(fitB.CLT_Burn$residuals)
  #fitBMC<-lm(B_sd_MCE/B_sd~B,data=stats1D)
  # fitB.CLT_Emp<-loess(B_ERR_CLT_Emp/B_RMSE~B,data=stats1D,span=0.5)
  # fitB.CLT_Emp
  # sd(fitB.CLT_Emp$residuals)
  #fitBLM<-lm(B_sd_LMI/B_sd~B,data=stats1D)
  fitB.LMI_Iter<-loess(B_ERR_LMI_Iter/B_RMSE~B,data=stats1D,span=0.5)
  fitB.LMI_Iter
  sd(fitB.LMI_Iter$residuals)
  fitB.LMI_Iter_Emp_2PL<-loess(B_ERR_LMI_Iter_Emp_2PL/B_RMSE~B,data=stats1D,span=0.5)
  fitB.LMI_Iter_Emp_2PL
  sd(fitB.LMI_Iter_Emp_2PL$residuals)
  #fitBo<-lm(B_sd_LMIo/B_sd~B,data=stats1D)
  fitB.LMI_Iter_Emp_2PNO<-loess(B_ERR_LMI_Iter_Emp_2PNO/B_RMSE~B,data=stats1D,span=0.5)
  fitB.LMI_Iter_Emp_2PNO
  sd(fitB.LMI_Iter_Emp_2PNO$residuals)
  #fitBMC<-lm(B_sd_MCE/B_sd~B,data=stats1D)
  fitB.Simple_2PL<-loess(B_ERR_Simple_2PL/B_RMSE~B,data=stats1D,span=0.5)
  fitB.Simple_2PL
  sd(fitB.Simple_2PL$residuals)
  #fitBLM<-lm(B_sd_LMI/B_sd~B,data=stats1D)
  fitB.Simple_2PNO<-loess(B_ERR_Simple_2PNO/B_RMSE~B,data=stats1D,span=0.5)
  fitB.Simple_2PNO
  sd(fitB.Simple_2PNO$residuals)
  # fitB.MCMC_Ann<-loess(B_ERR_MCMC_Ann/B_RMSE~B,data=stats1D,span=0.5)
  # fitB.MCMC_Ann
  # sd(fitB.MCMC_Ann$residuals)
  #fitBo<-lm(B_sd_LMIo/B_sd~B,data=stats1D)
  fitB.MCMC_Burn<-loess(B_ERR_MCMC_Burn/B_RMSE~B,data=stats1D,span=0.5)
  fitB.MCMC_Burn
  sd(fitB.MCMC_Burn$residuals)
  #fitBMC<-lm(B_sd_MCE/B_sd~B,data=stats1D)
  # fitB.MCMC_Emp<-loess(B_ERR_MCMC_Emp/B_RMSE~B,data=stats1D,span=0.5)
  # fitB.MCMC_Emp
  # sd(fitB.MCMC_Emp$residuals)
  ######plot2
  par(mar=c(4,6,4,m.b))
  seqB<-seq(range(stats1D$B)[1],range(stats1D$B)[2],length.out = 50)
  plot(stats1D$B,stats1D$B_ERR_MCMC_Burn/stats1D$B_RMSE,main=expression("Loess fit to"~hat(sigma)/sigma["RMSE"]),
       ylab="Ratio of Error Approximation to RMSE",xlab="Generated Intercept (2PNO)",ylim=YR,type="n")
  #lines(seqB,predict(fitB.CLT_Ann,newdata = data.frame(B=seqB)),lty=1)
  lines(seqB,predict(fitB.CLT_Burn,newdata = data.frame(B=seqB)),lty=2,col=2)
  #lines(seqB,predict(fitB.CLT_Emp,newdata = data.frame(B=seqB)),lty=3,col=3)
  #lines(seqB,predict(fitB.MCMC_Ann,newdata = data.frame(B=seqB)),lty=4,col=4,lwd=2)
  lines(seqB,predict(fitB.MCMC_Burn,newdata = data.frame(B=seqB)),lty=5,col=6,lwd=2)
  #lines(seqB,predict(fitB.MCMC_Emp,newdata = data.frame(B=seqB)),lty=6,col=1)
  legend("topright",c(#expression("CLT"["Ann"]),
    expression("CLT"["B"]),
    #expression("CLT"["Emp"]),
    #expression("MCMC"["Ann"]),
    expression("MCMC"["B"]) #,
    #expression("MCMC"["Emp"])
  ),#inset = 0.05,lty=1:6,col=c(1:4,6,1),lwd=c(1,1,1,2,2,1))
  inset = 0.05,lty=c(2,5),col=c(2,6),lwd=c(1,2))
  abline(h=1)
  par(mar=c(4+0.5,6,7,m.b))
  plot(stats1D$B,stats1D$B_ERR_MCMC_Burn/stats1D$B_RMSE,main=expression("Loess fit to Hessian"~hat(sigma)/sigma["RMSE"]),
       ylab="Ratio of Error Approximation to RMSE",xlab="Generated Intercept (2PNO)",ylim=YR,type="n")
  lines(seqB,predict(fitB.LMI_Iter,newdata = data.frame(B=seqB)),lty=1,col=1)
  lines(seqB,predict(fitB.LMI_Iter_Emp_2PL,newdata = data.frame(B=seqB)),lty=2,col=2)
  lines(seqB,predict(fitB.LMI_Iter_Emp_2PNO,newdata = data.frame(B=seqB)),lty=3,col=3,lwd=2)
  lines(seqB,predict(fitB.Simple_2PL,newdata = data.frame(B=seqB)),lty=4,col=4)
  lines(seqB,predict(fitB.Simple_2PNO,newdata = data.frame(B=seqB)),lty=5,col=6)
  legend("topright",c(expression("ICE"),
                      expression("IPCE"["L"]),
                      expression("IPCE"["O"]),
                      expression("SPCE"["L"]),
                      expression("SPCE"["O"])),inset = 0.05,lty=1:5,col=c(1:4,6),lwd=c(1,1,2,1,1))
  abline(h=1)
}