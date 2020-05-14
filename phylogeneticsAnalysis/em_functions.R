# Basic EM partitioning for DNA sequences

# Definition of functions:

# functions

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
    { # if a sequence has already been read, encountering a sequence header
      # means that this is the end of the previous sequence -> store it
      if(n_seq != 0)
      { seq_list[[n_seq]] <- seq_current }
      # a new sequence begins
      n_seq  <- n_seq + 1
      seq_current <- c()
    }
    else
    { seq_current <- c(seq_current, line) }
  }
  # store the last sequence
  seq_list[[n_seq]] <- seq_current
  # close file
  close(con)

  # store the sequences in a matrix
  sequences <- matrix(nrow=n_seq, ncol=length(seq_list[[1]]))
  for(i in 1:n_seq)
  { sequences[i,] <- seq_list[[i]] }

  return(sequences)
}

#' Converts a vector of characters containing a DNA sequence
#' into a vector of integers : A->1, C->2, G->3, T->4. Any
#' non ACGT character triggers an error.
#' @param sequence the DNA sequence stored as a vector of
#' characters.
#' @return a vector of integers.
#' @author Romain Groux
dna.to.int <- function(sequence)
{ seq.len <- length(sequence)
  seq.int <- vector(length=seq.len, mode="numeric")
  for(i in 1:seq.len)
      { if(sequence[i] == "A")
            { seq.int[i] <- 1 }
        else if(sequence[i] == "C")
            { seq.int[i] <- 2 }
        else if(sequence[i] == "G")
            { seq.int[i] <- 3 }
        else if(sequence[i] == "T")
            { seq.int[i] <- 4}
        else
            { stop(sprintf("Error! Unrecognized character in DNA sequence at position %d : %s", i, sequence[i])) }
    }
  return(seq.int)
}


expectation <- function(data, c, q) {

    N <- dim(data)[1]
    L <- dim(data)[2]
    K <- dim(c)[3]
    logc <- log(c+10**-100)
    logp <- matrix(0,nrow=N, ncol=K)
    for(b in 1:4) {
        x <- matrix(0,nrow=N, ncol=L)
        x[which(data == b, arr.ind=T)] <- 1
        logp <- logp + x %*% logc[b,,]
    }
    for(i in 1:N) {
        p[i,] <- q*exp(logp[i,] - max(logp[i,]))
        p[i,] <- p[i,] / sum(p[i,])
    }
    return(p)
}

maximization <- function(data, p, q) {

    N <- nrow(data)
    L <- ncol(data)
    K <- ncol(p)
    c <- array(data=0, dim=c(4, L, K))
    for(b in 1:4) {
        x <- matrix(0,nrow=N, ncol=L)
        x[which(data == b, arr.ind=T)] <- 1
        c[b,,] <- t(x) %*% p
    }
    for (i in 1:K) {c[,,i]=c[,,i] / (N*q[i])}
    return(c)
}

skew <- function(mat,f) {

    x <- mat**f
    mat <- t(t(x)/colSums(x))
    return(mat)
}


normalize <- function(c) {

    K <- dim(c)[3]
    for(i in 1:K) {
        comp <- rowMeans(c[,,i])
        ratio <- c[,,i]/comp
        c[,,i] <- t(t(ratio)/colSums(ratio))
    }
    return(c)
}

