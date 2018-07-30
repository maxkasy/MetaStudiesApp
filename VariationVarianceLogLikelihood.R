# tested both sym and asym versions
VariationVarianceLogLikelihood <-function(lambdabar, tauhat, betap, cutoffs, symmetric, X, sigma,normalizedsign) {
  ##show(X)
  n=length(X);
  betap=as.matrix(betap, length(betap),1);
  
  #%regressors for step function p
  TT=X/sigma;
  #disp(size(TT))
  
  Tpowers=matrix (0,n,length(cutoffs)+1);
  
  if (symmetric==1 ) {
    ##show(dim(Tpowers))
    ##show(length(TT))
    Tpowers[,1]=abs(TT)<cutoffs[1];
  
  
     if (length(cutoffs)>1) {
         for (m in (2:length(cutoffs))) {
              Tpowers[,m]=(abs(TT)<cutoffs[m])*(abs(TT)>=cutoffs[m-1]);
          }
          Tpowers[,length(cutoffs)+1]=abs(TT)>=cutoffs[length(cutoffs)];  
     }
    Tpowers[,length(cutoffs)+1]=abs(TT)>=cutoffs[length(cutoffs)];
    } else {
           
        Tpowers[,1]=TT<cutoffs[1];      
           
            
    if (length(cutoffs)>1) {
          for (m in 2:length(cutoffs)) {
              Tpowers[,m]=(TT<cutoffs[m])*(TT>=cutoffs[m-1]);  
          }
       
         
    }
        Tpowers[,length(cutoffs)+1]=(TT)>=cutoffs[length(cutoffs)]; 
  }
 
  

 
  
#  %%%%%%%%%%%%%%%%%%%%%%%%
#  % Calculate objective function LLH
  
#  %calculating components of the likelihood
#  %vector of estimated publication probabilities
 # ##show(Tpowers)
 # ##show(betap)
  phat=Tpowers%*%betap;
  
#@  %vector of un-truncated likelihoods
  if (normalizedsign==1) {
    fX=0.5*dnorm(X,lambdabar, sqrt(sigma^2 + tauhat^2))+0.5*dnorm(-X,lambdabar, sqrt(sigma^2 + tauhat^2));
  } else {
    fX=dnorm(X,lambdabar, sqrt(sigma^2 + tauhat^2));
  }
    
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
    ##show(normalizingconst)
    ##show('Phat')
    ##show(phat)
    ##show('fX')
    ##show(fX)
   # %cf equation 25 in the draft
    L<-phat*fX/normalizingconst;
    logL<-log(L);
    # %objective function; note the sign flip, since we are doing minimization
    LLH<--sum(log(L));
    ##show(logL)
    return(list(LLH=LLH,logL=logL))
    
  }
 
test<-FALSE
if (test) {
  cutoffs=c(1.96);
  betap = c(-1,0.099);
  symmetric<-1
  tauhat<-1
  normalizedsign<-1
  X<-1:10
  sigma<-1:10
  lambdabar<-1
  cdf<-VariationVarianceLogLikelihood(lambdabar, tauhat, betap, cutoffs, symmetric, X, sigma,normalizedsign)
}


