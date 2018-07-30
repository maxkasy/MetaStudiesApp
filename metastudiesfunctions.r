library(ggplot2)
# source("VariationVarianceLogLikelihood.R")
# source("RobustVariance.R")


metastudies_plot<-function(X,sigma,critval=1.96){
  
  n=length(X)
  significant<-(abs(X/sigma)>critval)
  nooutlier<-sigma<50;
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
      theme(aspect.ratio=1,
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "grey98", colour = NA))
  
}




metastudies_estimation=function(X,sigma, cutoffs,symmetric){
  #not sure I fully understand what's going on here...
  nn=length(X)
  if (symmetric == 0) {
    if (symmetric_p == 1) {
      LLH <- function (Psi) VariationVarianceLogLikelihood(Psi[1], Psi[2], c(1, Psi[-c(1,2)], rev(c(Psi[-c(1,2)])), 1),cutoffs,symmetric, X, sigma,0)
      Psihat0 = c(0,1, rep(0,length(cutoffs) + 1))
    } else {
      LLH <- function (Psi) VariationVarianceLogLikelihood(Psi[1], Psi[2], c(Psi[-c(1,2)],  1),cutoffs,symmetric, X, sigma,0)
      Psihat0 = c(0,1, 2*rep(0,length(cutoffs) + 1))
    }
  } else {
    LLH <-function (Psi)
      VariationVarianceLogLikelihood(0, Psi[1],c(Psi[-1], 1),cutoffs,symmetric, X, sigma,0)
  }
  LLH_only<-function (Psi){
    A<-LLH(Psi);
    return(A$LLH)
  }

  lower.b = c(-Inf,rep(0,length(Psihat0)-1))
  upper.b=rep(Inf,length(Psihat0))

  findmin<-nlminb(objective=LLH_only, start=Psihat0,lower=lower.b,upper=upper.b,control = list(eval.max = MaxEval, iter.max = MaxIter, abs.tol = Tol));
  Psihat<-findmin$par
  LLHmax<-findmin$objective

  Var_robust<-RobustVariance(stepsize, nn, Psihat, LLH,cluster_ID);
  se_robust<-sqrt(diag(Var_robust));

  #   SelectionTable(pathname,Psihat,se_robust,name,(identificationapproach==1),symmetric)

}


