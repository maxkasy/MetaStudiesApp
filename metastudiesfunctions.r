source("RobustVariance.r")


#parameters for optimization
MaxEval<-10^5
MaxIter<-10^5;
Tol<-10^(-8);
stepsize<-10^(-6);


Tpowers_fun=function(TT,cutoffs,symmetric){
  n=length(TT)
  Tpowers=matrix (0,n,length(cutoffs)+1)
  if (symmetric) TT=abs(TT)
  Tpowers[,1]=TT<cutoffs[1]
  if (length(cutoffs)>1) {
    for (m in 2:length(cutoffs)) {
      Tpowers[,m]=(TT<cutoffs[m])*(TT>=cutoffs[m-1]);  
    }
  }
  Tpowers[,length(cutoffs)+1]=(TT)>=cutoffs[length(cutoffs)]
  Tpowers    
}



VariationVarianceLogLikelihood <-function(lambdabar, tauhat, betap, cutoffs, symmetric, X, sigma, Tpowers) {
  n=length(X);
  betap=as.matrix(betap, length(betap),1);
  
  #  %vector of estimated publication probabilities
  phat=Tpowers%*%betap;
  
  #  %%%%%%%%%%%%%%%%%%%%%%%%
  #@  vector of un-truncated likelihoods
  fX=dnorm(X,lambdabar, sqrt(sigma^2 + tauhat^2));
  
  #  %normalizingconstant
  normalizingconst=rep(0,n);
  for (i in 1:n) {
    mu_i=lambdabar;
    sigma_tilde_i=sqrt((tauhat)^2 +sigma[i]^2);   
    prob_vec=rep(0,length(cutoffs)+1);
    if (symmetric==1) {
      for (m in 1:length(cutoffs)) {
        prob_vec[m+1]=(pnorm((cutoffs[m]*sigma[i]-mu_i)/sigma_tilde_i)-pnorm((-cutoffs[m]*sigma[i]-mu_i)/sigma_tilde_i)); 
      }
      prob_vec<-c(prob_vec,1);
      #prob_vec(end+1,1)=1;
      mean_Z1=prob_vec[-1]-prob_vec[-length(prob_vec)];
    } else {
      for (m in 1:length(cutoffs)) {
        prob_vec[m+1]=pnorm((cutoffs[m]*sigma[i]-mu_i)/sigma_tilde_i);
      }
      prob_vec<-c(prob_vec,1);
      mean_Z1<-prob_vec[-1]-prob_vec[-length(prob_vec)];
      
    }
    normalizingconst[i]<-t(mean_Z1)%*%betap;
    
  }
  
  #likelihood for each observation
  L<-phat*fX/normalizingconst;
  logL<-log(L);
  # objective function; note the sign flip, since we are doing minimization
  LLH<--sum(log(L));
  ##show(logL)
  return(list(LLH=LLH,logL=logL))
  
}







metastudies_estimation=function(X,sigma, cutoffs,symmetric){
    
  nn=length(X)
  #%regressors for step function p
  TT=X/sigma;
  Tpowers=Tpowers_fun(TT,cutoffs, symmetric)
  
  LLH <- function (Psi) VariationVarianceLogLikelihood(Psi[1], Psi[2], c(Psi[-c(1,2)],  1),cutoffs,symmetric, X, sigma, Tpowers)
  LLH_only<-function (Psi){
      A<-LLH(Psi);
      return(A$LLH)
  }  
  
  #starting values based on whether p is symmetric or not
  Psihat0 = c(0,1, rep(1,length(cutoffs)))



  lower.b = c(-Inf,rep(0,length(Psihat0)-1))
  upper.b=rep(Inf,length(Psihat0))

  
  
  findmin<-nlminb(objective=LLH_only, start=Psihat0,lower=lower.b,upper=upper.b,control = list(eval.max = MaxEval, iter.max = MaxIter, abs.tol = Tol));
  Psihat<-findmin$par
  LLHmax<-findmin$objective

  Var_robust<-RobustVariance(stepsize, nn, Psihat, LLH,1:nn);
  se_robust<-sqrt(diag(Var_robust));

  list(Psihat=Psihat, SE=se_robust)
}


