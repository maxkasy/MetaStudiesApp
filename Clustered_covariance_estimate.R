Clustered_covariance_estimate <-function(g,cluster_index) {
 # %given a matrix of moment condition values g, compute a clustering-robust
  #%estimate of the covariance matrix Sigma
  I<-order(cluster_index);
  cluster_index<-sort(cluster_index);
  g = g[I,]
  g = g - matrix(rep(apply(g,2,mean),length(I)),nrow=length(I),byrow=TRUE);
  gsum = apply(g,2,cumsum);
  index_diff <- cluster_index[-1] != cluster_index[-length(cluster_index)];
  index_diff <- c(index_diff,1);
  
  gsum = gsum[index_diff==1,]
  gsum=rbind(gsum[1,], diff(gsum));
  Sigma=1/(dim(g)[1]-1)*(t(gsum)%*%gsum);
                         
  
  return (Sigma)
}