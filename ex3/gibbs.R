library(spatial)
library(MASS)

o <- as.matrix(read.table('seismic.dat'))
o <- matrix(o, 75, 75)

drawL <- function(o) {
  Lsand  <- dnorm(o, mean=0.02, sd=0.06)
  Lshale <- dnorm(o, mean=0.08, sd=0.06)
  Pshale <- Lshale / (Lsand + Lshale)
  r <- matrix(runif(o), nrow(o), ncol(o))
  return(Pshale > r)
}

complit <- as.matrix(read.table('complit.dat'))

modF <- function(i, n) {
  # i is from 0 to n+1
  # should return from 1 to n
  i <- i - 1
  # i is now from -1 to n
  i <- i %% n
  # i is now from 0 to n-1
  i <- i + 1
  # i is now from 1 to n
  return(i)
}

sumAround <- function(L, r, c, comp) {
  s <- 0
  s <- s + (L[modF(r+1,nrow(L)),c] == comp)
  s <- s + (L[modF(r-1,nrow(L)),c] == comp)
  s <- s + (L[r,modF(c+1,ncol(L))] == comp)
  s <- s + (L[r,modF(c-1,ncol(L))] == comp)
}

pseudoProb <- function(L, beta) {
  tot <- 0
  
  for (i in 1:nrow(L)) {
    for (j in 1:ncol(L)) {
      tot <- tot + beta * sumAround(L, i, j, L[i,j]) 
      tot <- tot - log(exp(beta * sumAround(L, i, j, 0))
                       + exp(beta * sumAround(L, i, j, 1)))
    }
  }
  
  return(tot)
}

posterior <- function(L, obs, beta) {
  tot <- 0
  for (i in 1:nrow(obs)) {
    for (j in 1:ncol(obs)) {
      if (L[i,j] == 0) {
        tot <- tot - (obs[i,j] - 0.02)^2 / (2*0.06^2)
      } else {
        tot <- tot - (obs[i,j] - 0.08)^2 / (2*0.06^2)
      }
      
      tot <- tot + beta * (L[i,j] == L[i,modF(j+1, ncol(obs))])
      tot <- tot + beta * (L[i,j] == L[modF(i+1, nrow(obs)),j])
    }
  }
  return(tot)
}

gibbs <- function(obs, beta, iter) {
  
  # Resulting images stacked behind eachother
  res <- array(NA, c(nrow(obs), ncol(obs), iter))
  
  # Initial value, something else than just 0, to get it going
  res[,,1] = 1 * (obs > 0)
  
  p0p <- dnorm(obs, mean=0.02, sd=0.06)
  p1p <- dnorm(obs, mean=0.08, sd=0.06)
  
  for (it in 2:iter) {
    cat("Iter: ", it, "\n")
    # Image of most current values to use in Gibbs update
    use <- res[,,it-1]
    
    # Could (should?) iterate randomly here for example
    for (i in 1:nrow(obs)) {
      for (j in 1:ncol(obs)) {
        
        # Proportional to the probability of 0 and 1 in this square
        sum0 <- sumAround(use, i, j, 0)
        p0q <- p0p[i,j] * exp(beta * sum0)
        p1q <- p1p[i,j] * exp(beta * (4 - sum0))
        
        # Probability of 0 in this square
        p0 <- p0q / (p0q + p1q)
        
        # Insert 0 or 1?
        r <- runif(1)
        if (r < p0) {
          res[i,j,it] = 0
          use[i,j] = 0
        } else {
          res[i,j,it] = 1
          use[i,j] = 1
        }
      }
    }

    if (it %% 10 == 0) {
      image(res[,,it])
    }
  }
  return(res)
}