
par(mfrow=c(2,2))

pnts1<-list(y=c(0.04+0.025*1:5+rnorm(5,0.05,0.025),0.7+0.04*1:5+rnorm(5,0.05,0.03)),
           x=c(0.05+0.1*0:9))
pnts1$theta=qnorm(pnts1$x)

plot(pnts1$x,pnts1$y,xlim=c(0,1),ylim=c(0,1),xlab="score deciles",ylab="P(Y=1|score decile)",main="High Discrimination, Mean Difficulty")
lines(pnts1$x,pnts1$y)


pnts2<-list(y=c(0.45+0.04*1:10+rnorm(10,0.05,0.02)),
            x=c(0.05+0.1*0:9))
pnts2$theta=qnorm(pnts2$x)

plot(pnts2$x,pnts2$y,xlim=c(0,1),ylim=c(0,1),xlab="score deciles",ylab="P(Y=1|score decile)",main="Low Discrimination, Low Difficulty")
lines(pnts2$x,pnts2$y)


pnts3<-list(y=c(0.05+0.02*1:8+rnorm(8,0.05,0.025),0.8+0.04*1:2+rnorm(2,0.05,0.03)),
            x=c(0.05+0.1*0:9))
pnts3$theta=qnorm(pnts3$x)

plot(pnts3$x,pnts3$y,xlim=c(0,1),ylim=c(0,1),xlab="score deciles",ylab="P(Y=1|score decile)",main="High Discrimination, High Difficulty")
lines(pnts3$x,pnts3$y)

pnts4<-list(y=c(0.45+0.02*1:10+rnorm(10,0.02,0.06)),
            x=c(0.05+0.1*0:9))
pnts4$theta=qnorm(pnts4$x)

plot(pnts4$x,pnts4$y,xlim=c(0,1),ylim=c(0,1),xlab="score deciles",ylab="P(Y=1|score decile)",main="Uninformative")
lines(pnts4$x,pnts4$y)


irt3fit<-function(thet,gam=0.1,alph=1,bet=0) {
  gam+(1-gam)*1/(1+exp(-1.7*alph*(thet-bet))) 
}
ogive3fit<-function(thet,gam=0.1,alph=1,bet=0) {
  gam+(1-gam)*pnorm(alph*thet-bet)
}

theta<-qnorm(0.005+0:99*0.01)

fit1<-nls(formula=y~irt3fit(thet,gam,alph,bet),data=data.frame(y=pnts1$y,thet=pnts1$theta),start = list(gam=0.1,alph=1,bet=0))

plot(pnts1$theta,pnts1$y,xlim=c(-3,3),ylim=c(0,1),xlab="score deciles as Z-scores",ylab="P(Y=1|Z-score)",main="IRT-3PL: High Discrimination, Mean Difficulty")
lines(pnts1$theta,pnts1$y)
lines(theta,predict(fit1,newdata = data.frame(thet=theta)),lwd=2,col=2,lty=3)
eval(parse(text=paste0('text(2,0.55,expression(paste(hat(alpha)," = ",',format(fit1$m$getPars()[2],digits=2),')))')))
eval(parse(text=paste0('text(2,0.4,expression(paste(hat(beta)," = ",',format(fit1$m$getPars()[3],digits=2),')))')))
eval(parse(text=paste0('text(2,0.25,expression(paste(hat(gamma)," = ",',format(fit1$m$getPars()[1],digits=2),')))')))

fit2<-nls(formula=y~irt3fit(thet,gam,alph,bet),data=data.frame(y=pnts2$y,thet=pnts2$theta),start = list(gam=0.1,alph=1,bet=0))

plot(pnts2$theta,pnts2$y,xlim=c(-3,3),ylim=c(0,1),xlab="score deciles as Z-scores",ylab="P(Y=1|Z-score)",main="IRT-3PL: Low Discrimination, Low Difficulty")
lines(pnts2$theta,pnts2$y)
lines(theta,predict(fit2,newdata = data.frame(thet=theta)),lwd=2,col=2,lty=3)
eval(parse(text=paste0('text(2,0.55,expression(paste(hat(alpha)," = ",',format(fit2$m$getPars()[2],digits=2),')))')))
eval(parse(text=paste0('text(2,0.4,expression(paste(hat(beta)," = ",',format(fit2$m$getPars()[3],digits=2),')))')))
eval(parse(text=paste0('text(2,0.25,expression(paste(hat(gamma)," = ",',format(fit2$m$getPars()[1],digits=2),')))')))

fit3<-nls(formula=y~irt3fit(thet,gam,alph,bet),data=data.frame(y=pnts3$y,thet=pnts3$theta),start = list(gam=0.1,alph=1,bet=0))

plot(pnts3$theta,pnts3$y,xlim=c(-3,3),ylim=c(0,1),xlab="score deciles as Z-scores",ylab="P(Y=1|Z-score)",main="IRT-3PL: High Discrimination, High Difficulty")
lines(pnts3$theta,pnts3$y)
lines(theta,predict(fit3,newdata = data.frame(thet=theta)),lwd=2,col=2,lty=3)
eval(parse(text=paste0('text(2,0.55,expression(paste(hat(alpha)," = ",',format(fit3$m$getPars()[2],digits=2),')))')))
eval(parse(text=paste0('text(2,0.4,expression(paste(hat(beta)," = ",',format(fit3$m$getPars()[3],digits=2),')))')))
eval(parse(text=paste0('text(2,0.25,expression(paste(hat(gamma)," = ",',format(fit3$m$getPars()[1],digits=2),')))')))

fit4<-nls(formula=y~irt3fit(thet,gam,alph,bet),data=data.frame(y=pnts4$y,thet=pnts4$theta),start = list(gam=0.1,alph=1,bet=0))

plot(pnts4$theta,pnts4$y,xlim=c(-3,3),ylim=c(0,1),xlab="score deciles as Z-scores",ylab="P(Y=1|Z-score)",main="IRT-3PL: Uninformative")
lines(pnts4$theta,pnts4$y)
lines(theta,predict(fit4,newdata = data.frame(thet=theta)),lwd=2,col=2,lty=3)
eval(parse(text=paste0('text(2,0.55,expression(paste(hat(alpha)," = ",',format(fit4$m$getPars()[2],digits=2),')))')))
eval(parse(text=paste0('text(2,0.4,expression(paste(hat(beta)," = ",',format(fit4$m$getPars()[3],digits=2),')))')))
eval(parse(text=paste0('text(2,0.25,expression(paste(hat(gamma)," = ",',format(fit4$m$getPars()[1],digits=2),')))')))

############ Ogive fits

fit1<-nls(formula=y~ogive3fit(thet,gam,alph,bet),data=data.frame(y=pnts1$y,thet=pnts1$theta),start = list(gam=0.0,alph=2,bet=0))

plot(pnts1$theta,pnts1$y,xlim=c(-3,3),ylim=c(0,1),xlab="score deciles as Z-scores",ylab="P(Y=1|Z-score)",main="3PNO: High Discrimination, Mean Difficulty")
lines(pnts1$theta,pnts1$y)
lines(theta,predict(fit1,newdata = data.frame(thet=theta)),lwd=2,col=2,lty=3)
eval(parse(text=paste0('text(2,0.55,expression(paste(hat(a)," = ",',format(fit1$m$getPars()[2],digits=2),')))')))
eval(parse(text=paste0('text(2,0.4,expression(paste(hat(b)," = ",',format(fit1$m$getPars()[3],digits=2),')))')))
eval(parse(text=paste0('text(2,0.25,expression(paste(hat(c)," = ",',format(fit1$m$getPars()[1],digits=2),')))')))

fit2<-nls(formula=y~ogive3fit(thet,gam,alph,bet),data=data.frame(y=pnts2$y,thet=pnts2$theta),start = list(gam=0.1,alph=1,bet=0))

plot(pnts2$theta,pnts2$y,xlim=c(-3,3),ylim=c(0,1),xlab="score deciles as Z-scores",ylab="P(Y=1|Z-score)",main="3PNO: Low Discrimination, Low Difficulty")
lines(pnts2$theta,pnts2$y)
lines(theta,predict(fit2,newdata = data.frame(thet=theta)),lwd=2,col=2,lty=3)
eval(parse(text=paste0('text(2,0.55,expression(paste(hat(a)," = ",',format(fit2$m$getPars()[2],digits=2),')))')))
eval(parse(text=paste0('text(2,0.4,expression(paste(hat(b)," = ",',format(fit2$m$getPars()[3],digits=2),')))')))
eval(parse(text=paste0('text(2,0.25,expression(paste(hat(c)," = ",',format(fit2$m$getPars()[1],digits=2),')))')))

fit3<-nls(formula=y~ogive3fit(thet,gam,alph,bet),data=data.frame(y=pnts3$y,thet=pnts3$theta),start = list(gam=0.1,alph=1,bet=0))

plot(pnts3$theta,pnts3$y,xlim=c(-3,3),ylim=c(0,1),xlab="score deciles as Z-scores",ylab="P(Y=1|Z-score)",main="3PNO: High Discrimination, High Difficulty")
lines(pnts3$theta,pnts3$y)
lines(theta,predict(fit3,newdata = data.frame(thet=theta)),lwd=2,col=2,lty=3)
eval(parse(text=paste0('text(2,0.55,expression(paste(hat(a)," = ",',format(fit3$m$getPars()[2],digits=2),')))')))
eval(parse(text=paste0('text(2,0.4,expression(paste(hat(b)," = ",',format(fit3$m$getPars()[3],digits=2),')))')))
eval(parse(text=paste0('text(2,0.25,expression(paste(hat(c)," = ",',format(fit3$m$getPars()[1],digits=2),')))')))

fit4<-nls(formula=y~ogive3fit(thet,gam,alph,bet),data=data.frame(y=pnts4$y,thet=pnts4$theta),start = list(gam=0.5,alph=0.2,bet=0))

plot(pnts4$theta,pnts4$y,xlim=c(-3,3),ylim=c(0,1),xlab="score deciles as Z-scores",ylab="P(Y=1|Z-score)",main="3PNO: Uninformative")
lines(pnts4$theta,pnts4$y)
lines(theta,predict(fit4,newdata = data.frame(thet=theta)),lwd=2,col=2,lty=3)
eval(parse(text=paste0('text(2,0.55,expression(paste(hat(a)," = ",',format(fit4$m$getPars()[2],digits=2),')))')))
eval(parse(text=paste0('text(2,0.4,expression(paste(hat(b)," = ",',format(fit4$m$getPars()[3],digits=2),')))')))
eval(parse(text=paste0('text(2,0.25,expression(paste(hat(c)," = ",',format(fit4$m$getPars()[1],digits=2),')))')))
