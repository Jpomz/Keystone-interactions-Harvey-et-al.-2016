jacobian_with_cannibals <- function(m){
  # m is the adjacency matrix of the graph, dim(m) = S*S
  # undirected interactions (or bidirectional) are such that m_ij = m_ji = 1
  # directected interactions (such as antagonistic ones) are such that m_ij = 1 and m_ji = 0 if i feeds on j (j -> i)
  # S being the total number of species
  
  if (dim(m)[1] == dim(m)[2]){ # Is m a square matrix?
    J <- m
    J[which(J < t(J))] <- -t(J)[which(J < t(J))]
    
    L <- (sum(abs(J)) - sum(abs(diag(J))))/2 # here the links doesn't comprise the cannibalism (remove the links on the diagonal)
    strength <- rnorm(L, sd = 0.1) + 1 # L values for interaction strength magnitude drawn from a normal distribution
    upper_tri_str_index <- which((upper.tri(J) == TRUE) & (J != 0), arr.ind = TRUE) # non-zero elements in the upper triangle matrix
    
    lower_tri_str_index <- cbind(upper_tri_str_index[, 2], upper_tri_str_index[, 1]) # so indices match between upper triangle and lower triangle elements
    J[upper_tri_str_index] <- J[upper_tri_str_index]*strength
    
    strength <- rnorm(L, sd = 0.1) + 1
    J[lower_tri_str_index] <- J[lower_tri_str_index]*strength
    
    # parameterize the diagonal (negative interaction strengths due to cannibalism)
    diag(J) <- diag(J) * (rnorm(length(diag(J)), sd = 0.1) + 1)* (-1)
    return(J)
  }
  else {
    stop("m is not a square matrix.")
  }  
}