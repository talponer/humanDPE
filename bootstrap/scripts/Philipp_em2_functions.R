# probabilistic partitioning with EM - extended algorithm
# Definition of functions:

  expectation = function(data, c, q) {

    N=dim(data)[1]; L=dim(data)[2]; K=dim(c)[3]
    logc=log(c+10**-100)
    logp=matrix(0,nrow=N, ncol=K)
    for(b in 1:4) {
      x = matrix(0,nrow=N, ncol=L)
      x[which(data == b, arr.ind=T)] =1
      logp = logp + x %*% logc[b,,]
      }
    p=matrix(0,nrow=N, ncol=K)
    for(i in 1:N) {
      p[i,] = q*exp(logp[i,] - max(logp[i,]))
      p[i,] = p[i,] / sum(p[i,])
      }
    return(p)
    }

  expectation2= function(data, c, q) {

    N=dim(data)[1]; L=dim(data)[2]; K=dim(c)[3]

    DATA=matrix(nrow=N,ncol=0)
    for(j in 1:L) {for(b in 1:4) {DATA=cbind(DATA,as.numeric(data[,j] == b))}}
       
    logc=log(c+10**-100)   
    logp=matrix(0,nrow=N, ncol=K) 

    for(i in 1:N) {for(k in 1:K) {logp[i,k]=sum(DATA[i,]*as.vector(logc[,,k]))}} 

    p=matrix(0,nrow=N, ncol=K) 
    for(i in 1:N) {
      p[i,] = q*exp(logp[i,] - max(logp[i,]))
      p[i,] = p[i,] / sum(p[i,])
      }
    return(p)
    }

  maximization = function(data, p, q) {

    N = nrow(data); L=ncol(data); K=ncol(p)
    c = array(data=0, dim=c(4, L, K)) 
    for(b in 1:4) {
      x = matrix(0,nrow=N, ncol=L)
      x[which(data == b, arr.ind=T)] =1
      c[b,,] = t(x) %*% p 
      }
    for (i in 1:K) {c[,,i]=c[,,i] / (N*q[i])}
    return(c)
    }

  skew = function(mat,f) {
     x = mat**f
     mat = t(t(x)/colSums(x))
     return(mat)
  }

  normalize = function(c, w, comp) { 

    K=dim(c)[3]
    for(i in 1:K) {
      comp=rowMeans(c[,,i])
      ratio=c[,,i]/comp
      c[,,i]=t(t(ratio+w)/colSums(ratio+w))
      }
    return(c)
    }

