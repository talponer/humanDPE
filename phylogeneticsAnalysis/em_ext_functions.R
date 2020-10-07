#' Reads a FASTA file and returns a character matrix
#' with the sequences on the rows.
#' @param file_fasta the path to file to read.
#' @return a matrix of characters where each row contains
#' a sequence from the file.
#' @author Romain Groux

read.fasta <- function(file_fasta)
{
    con      <- file(file_fasta, "r")
    reading  <- TRUE
    n_seq    <- 0
    seq_list <- list()
    while(reading)
    {
        line  <- strsplit(readLines(con, n=1), split='')
        
        if(length(line) == 0)
        { reading <- FALSE
            next
        }
        
        line <- line[[1]]
        
        if(line[1] == '>')
        { ## if a sequence has already been read, encountering a
          ## sequence header means that this is the end of the previous
          ## sequence -> store it
            if(n_seq != 0)
            { seq_list[[n_seq]] <- seq_current }
            ## a new sequence begins
            n_seq  <- n_seq + 1
            seq_current <- c()
        }
    else
        { seq_current <- c(seq_current, line) }
    }
    ## store the last sequence
    seq_list[[n_seq]] <- seq_current
    ## close file
    close(con)

    ## store the sequences in a matrix
    sequences <- matrix(nrow=n_seq, ncol=length(seq_list[[1]]))
    for(i in 1:n_seq)
    { sequences[i,] <- seq_list[[i]] }
    
    return(sequences)
}


## probabilistic partitioning with EM - extended algorithm
## Definition of functions:

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

