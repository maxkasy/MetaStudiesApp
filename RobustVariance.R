# tested
source("Clustered_covariance_estimate.R")
RobustVariance<-function(stepsize,n, thetahat, LLH,cluster_ID) {
  Info <- matrix(0,length(thetahat),length(thetahat));
  for (n1 in 1:length(thetahat)) {
    for (n2 in 1:length(thetahat)) {
      thetaplusplus<-thetahat;
      thetaplusminus<-thetahat;
      thetaminusplus<-thetahat;
      thetaminusminus<-thetahat;
      
      
      
      thetaplusplus[n1]<-thetaplusplus[n1]+stepsize;
      thetaplusplus[n2]<-thetaplusplus[n2]+stepsize;
      LLH_plusplus<-LLH(thetaplusplus);
      LLH_plusplus<-LLH_plusplus$LLH
      
      thetaplusminus[n1]<-thetaplusminus[n1]+stepsize;
      thetaplusminus[n2]<-thetaplusminus[n2]-stepsize;
      LLH_plusminus<-LLH(thetaplusminus);
      LLH_plusminus<-LLH_plusminus$LLH
      
      thetaminusplus[n1]<-thetaminusplus[n1]-stepsize;
      thetaminusplus[n2]<-thetaminusplus[n2]+stepsize;
      LLH_minusplus<-LLH(thetaminusplus);
      LLH_minusplus<-LLH_minusplus$LLH
      
      thetaminusminus[n1]<-thetaminusminus[n1]-stepsize;
      thetaminusminus[n2]<-thetaminusminus[n2]-stepsize;
      LLH_minusminus<-LLH(thetaminusminus);
      LLH_minusminus<-LLH_minusminus$LLH
      
      Info[n1,n2]=((LLH_plusplus-LLH_plusminus)/(2*stepsize)-(LLH_minusplus-LLH_minusminus)/(2*stepsize))/(2*stepsize);
    }
  }
  
  Var=solve(Info);
  #show(Var)
  #%Calculate misspecification-robust standard errors
  score_mat <- matrix(0,n,length(thetahat));
  for (n1 in 1:length(thetahat)) {
    theta_plus<-thetahat;
    theta_plus[n1]<-theta_plus[n1]+stepsize;
    funvalue<-LLH(theta_plus);
    
    #funvalue<-funvalue$LLH;
    LLH_plus<-funvalue$LLH;
    logL_plus<-funvalue$logL;
    
    theta_plus<-thetahat;
    theta_plus[n1]<-theta_plus[n1]-stepsize;
    funvalue<-LLH(theta_plus);
   # funvalue<-funvalue$LLH;
    LLH_minus<-funvalue$LLH;
    logL_minus<-funvalue$logL;
    
    score_mat[,n1]=(logL_plus-logL_minus)/(2*stepsize);
  }
  Cov=Clustered_covariance_estimate(score_mat,cluster_ID);
  #show(score_mat)
  #show(Cov)
  Var_robust=n*solve(Info)%*%Cov%*%solve(Info);
  
return(Var_robust)
}