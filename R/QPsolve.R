# Performs whole genome breed composition prediction.
# 
# @param Y numeric vector of genotypes (with names as SNPs) from a single animal.
#   coded as dosage of allele B {0, 1, 2}
# @param X numeric matrix of allele frequencies from reference animals
# @param p numeric indicating number of breeds represented in X
# @param names character names of breeds
# @return data.frame of breed composition estimates
# @import quadprog
# @export
QPsolve <- function(Y, X, p = 4, names = c("Duroc", "Hampshire", "Landrace", "Yorkshire")) {
  
  # Remove NAs from Y and remove corresponding
  #   SNPs from X. Ensure Y is numeric
  Ymod <- Y[!is.na(Y)]
  Xmod <- X[names(Ymod), ]
  
  # perfom steps needed to solve OLS by framing
  # as a QP problem
  # Rinv - matrix should be of dimensions px(p+1) where p is the number of variables in X
  Rinv <- solve(chol(t(Xmod) %*% Xmod))
  
  # C - the first column is a sum restriction (all equal to 1) and the rest of the columns an identity matrix
  C <- cbind(rep(1, p), diag(p))
  
  # b2 - This should be a vector of length p+1 the first element is the value of the sum (1)
  #   the other elements are the restriction of individual coefficients (>)
  #   so a value 0 produces positive coefficients
  b2 <- c(1, rep(0, p))
  
  # dd - this should be a matrix NOT a vector
  dd <- (t(Ymod) %*% Xmod)
  
  qp <- solve.QP(Dmat = Rinv, factorized = TRUE, dvec = dd, Amat = C, bvec = b2, meq = 1)
  beta <- qp$solution
  rr <- cor(Ymod, Xmod %*% beta)^2
  result <- c(beta, rr)
  names(result) <- c(names, "R2")
  return(result)
}