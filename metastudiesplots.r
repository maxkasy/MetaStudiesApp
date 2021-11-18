library(ggplot2)
library(reshape2)

#for drawing
critval=1.96

metastudies_plot<-function(X,sigma){
  n=length(X)
  significant<-(abs(X/sigma)>critval)
  nooutlier= (sigma<30*mean(sigma))&(abs(X) < 30*mean(abs(X)));
  dat<-data.frame(X, sigma, as.factor(significant&nooutlier))
  names(dat)=c("xvar", "yvar","significant")
  rangeX=1.1*max(max(abs(X)), max(abs(sigma[nooutlier]))*critval)
  
  dat<-dat[order(dat$significant),]
  
  ggplot(dat, aes(x=xvar,y=yvar)) +
    xlab("X")+
    ylab(expression(sigma))+
    geom_abline(intercept = 0,slope=1/critval,color="grey")+ 
    geom_abline(intercept = 0,slope=-1/critval,color="grey")+
    geom_point(size = 4,aes(colour = significant,
                            fill = significant), alpha=min(.8,max(40/n,.3)))+ 
    #scale_fill_manual(values=c("grey", "blue")) + 
    scale_colour_manual(values=c("grey50", "blue")) +
    scale_x_continuous(expand = c(0,0),limits = c(-rangeX,rangeX))+
    scale_y_continuous(expand = c(0,0), limits = c(0,rangeX/critval))+
    theme(legend.position="top",
          #aspect.ratio=1,
          #panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "grey95", colour = NA))
  
}

z_histogram=function(X,sigma){
  Z=X/sigma
  n=length(Z)
  ll=floor(min(Z));
  uu=ceiling(max(Z));
  
    if (n>=30) {
      uu2<-ceiling(max((uu-.36)/.32,0))*.32+.36;
      ll2<-floor(min((ll+.36)/.32,0))*.32-.36;
      edges<-c(seq(from=ll2,
                   to=-0.36,
                   by=0.32), 0, seq(from=0.36,
                                    to=uu2,
                                    by=0.32));
    } else {
      uu2<-ceiling(max((uu-.68)/.64,0))*.64+.68;
      ll2<-floor(min((ll+.68)/.64,0))*.64-.68;
      edges<-c(seq(from=ll2,
                   to=-0.68,
                   by=0.64), 0, seq(from=0.68,
                                    to=uu2,
                                    by=0.64));
    }

  
  ggplot(data = as.data.frame(Z), aes(Z))+
    geom_histogram(aes(y = ..density..),
                   fill = 'blue',
                   breaks=edges)+
    geom_vline(xintercept =-1.96,color='grey')+
    geom_vline(xintercept =1.96, color='grey')+
    xlab('Z')+
    ylab('Density')+
    xlim(c(min(edges),max(edges)))+
    theme(#panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "grey95", colour = NA))
  
}


estimates_plot<-function(X, sigma, cutoffs, symmetric, estimates, model="normal"){
  
  n=500
  Psihat=estimates$Psihat
  rangeZ=3 
  dens=data.frame(z=seq(-rangeZ,rangeZ,length.out =n))
  shift=as.integer(model=="t")
  
  Tpowers=Tpowers_fun(dens$z,cutoffs,symmetric)
  betap=as.vector(c(Psihat[-(1:(2+shift))],  1))
  dens$p=Tpowers%*%betap  
  
  if (model=="t") df=Psihat[3]
    else df=Inf
  
  dens$f=dt(((dens$z - Psihat[1])/ Psihat[2]), df=df)/Psihat[2]
  names(dens)[names(dens) == 'f'] <- 'density of true effect'
  names(dens)[names(dens) == 'p'] <- 'publication probability'
  
  dens=melt(dens, id="z")
  ggplot(dens, aes(x=z, y=value)) +
    xlab(paste("Z, ", intToUtf8(952)))+
    geom_line(size=2, color="blue") +
    facet_grid(variable ~ .,  scales = "free_y") +
    expand_limits(y = 0) +
    scale_x_continuous(breaks =-3:3) +
    theme(#panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "grey95", colour = NA))
  

}