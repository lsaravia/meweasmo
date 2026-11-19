#' R implementation of generateGLVparmsFromAdj
#' 
#' Generates random parameters of a generalized Lotka-Volterra model from an adjacency matrix. 
#' It is equivalent to [generateGLVparmsFromAdj()]: it has the same arguments and, for the same 
#' random seed, it gives the same results (this is checked by the package tests). 
#' It is written in R, so it is slower but easier to read and modify, and it can be used as a reference 
#' implementation. See [generateGLVparmsFromAdj()] for a detailed description of how the 
#' parameters are generated.
#' 
#' @inheritParams generateGLVparmsFromAdj
#' @inherit generateGLVparmsFromAdj return
#' @seealso [generateGLVparmsFromAdj()] for the Rcpp implementation, [generateRandomGLVadjMat()] to 
#'          generate `adjM`, and [metaWebNetAssemblyGLV()].
#' @examples
#' set.seed(1)
#' adj <- generateRandomGLVadjMat(20, 0.2, c(0.1, 0.3, 0.3, 0.15, 0.15))
#' set.seed(2)
#' a <- generateGLVparmsFromAdj(adj, ef = 0.1)
#' set.seed(2)
#' b <- generateGLVparmsFromAdjR(adj, ef = 0.1)
#' all.equal(a, b)
#' @importFrom stats, rbinom, rgamma, runif
#' @export
generateGLVparmsFromAdjR <- function(adjM, ef, predIntAvg = 0.01, selfLimAvg = c(0.01, 0.01, 0.01),
                                    migrAvg = 0.0, preserveInt = 0, predIntSd = 0.01, rndType = 0) {
  rho <- nrow(adjM)  # metacommunity size
  
  A <- matrix(FALSE,rho,rho)  # Local adjacency matrix
  r <- numeric(rho)  # Growth rates
  m <- numeric(rho)  # Migration rates
  Bas <- rep(TRUE, rho)  # Basal species
  Mut <- rep(FALSE, rho)  # Mutualistic species
  
  if (rndType == 0) {
    # Uniform distribution
    a_predInt <- 0
    b_predInt <- predIntAvg * 2
  } else if (rndType == 1) {
    # Gamma distribution
    shape <- (predIntAvg / predIntSd)^2
    rate <- predIntAvg / predIntSd^2  # rate = 1/scale, scale = sd^2/mean
  }
  
  # Fill adjacency matrix A and number of species SL
  diag(A) <- TRUE
  if (preserveInt == 0) {
    for (i in 1:rho) {
      for (j in 1:rho) {
        if (adjM[i, j] != 0 && !A[i, j]) {
          if (adjM[i, j] < 0 && adjM[j, i] > 0) {  # j is the predator, i is the prey
            A[i, j] <- A[j, i] <- TRUE
            if (rndType == 0) {
              adjM[i, j] <- -runif(1, a_predInt, b_predInt)
            } else if (rndType == 1) {
              adjM[i, j] <- -rgamma(1, shape, rate)
            }
            adjM[j, i] <- ef * -adjM[i, j]
            Bas[j] <- FALSE  # Predators are not basal
          } else if (adjM[i, j] < 0 && adjM[j, i] < 0) {  # Competition
            A[i, j] <- A[j, i] <- TRUE
            if (rndType == 0) {
              adjM[i, j] <- -runif(1, a_predInt, b_predInt)
              adjM[j, i] <- -runif(1, a_predInt, b_predInt)
            } else if (rndType == 1) {
              adjM[i, j] <- -rgamma(1, shape, rate)
              adjM[j, i] <- -rgamma(1, shape, rate)
            }
          } else if (adjM[i, j] > 0 && adjM[j, i] > 0) {  # Mutualism
            A[i, j] <- A[j, i] <- TRUE
            if (rndType == 0) {
              adjM[i, j] <- runif(1, a_predInt, b_predInt)
              adjM[j, i] <- runif(1, a_predInt, b_predInt)
            } else if (rndType == 1) {
              adjM[i, j] <- rgamma(1, shape, rate)
              adjM[j, i] <- rgamma(1, shape, rate)
            }
            Mut[j] <- Mut[i] <- TRUE
          } else if (adjM[i, j] > 0 && adjM[j, i] == 0) {  # Commensalism
            A[i, j] <- TRUE
            if (rndType == 0) {
              adjM[i, j] <- runif(1, a_predInt, b_predInt)
            } else if (rndType == 1) {
              adjM[i, j] <- rgamma(1, shape, rate)
            }
            Bas[j] <- FALSE
          } else if (adjM[i, j] < 0 && adjM[j, i] == 0) {  # Amensalism
            A[i, j] <- TRUE
            if (rndType == 0) {
              adjM[i, j] <- -runif(1, a_predInt, b_predInt)
            } else if (rndType == 1) {
              adjM[i, j] <- -rgamma(1, shape, rate)
            }
          }
        }
      }
    }
  } else if (preserveInt == 1) {  # Do not touch the interaction values
    for (i in 1:rho) {
      for (j in 1:rho) {
        if (adjM[i, j] != 0 && !A[i, j]) {
          if (adjM[i, j] < 0 && adjM[j, i] > 0) {  # j is the predator, i is the prey
            A[i, j] <- A[j, i] <- TRUE
            Bas[j] <- FALSE  # Predators are not basal
          } else if (adjM[i, j] < 0 && adjM[j, i] < 0) {  # Competition
            A[i, j] <- A[j, i] <- TRUE
          } else if (adjM[i, j] > 0 && adjM[j, i] > 0) {  # Mutualism
            A[i, j] <- A[j, i] <- TRUE
            Mut[j] <- Mut[i] <- TRUE
          } else if (adjM[i, j] > 0 && adjM[j, i] == 0) {  # Commensalism
            A[i, j] <- TRUE
            Bas[j] <- FALSE
          } else if (adjM[i, j] < 0 && adjM[j, i] == 0) {  # Amensalism
            A[i, j] <- TRUE
          }
        }
      }
    }
  } else if (preserveInt == 2) {
    for (i in 1:rho) {
      for (j in 1:rho) {
        if (adjM[i, j] != 0 && !A[i, j]) {
          if (adjM[i, j] < 0 && adjM[j, i] > 0) {  # j is the predator, i is the prey
            A[i, j] <- A[j, i] <- TRUE
            adjM[i, j] <- -runif(1, 0, -2 * adjM[i, j])
            adjM[j, i] <- ef * -adjM[i, j]
            Bas[j] <- FALSE  # Predators are not basal
          } else if (adjM[i, j] < 0 && adjM[j, i] < 0) {  # Competition
            A[i, j] <- A[j, i] <- TRUE
            adjM[i, j] <- -runif(1, 0, -2 * adjM[i, j])
            adjM[j, i] <- -runif(1, 0, -2 * adjM[j, i])
          } else if (adjM[i, j] > 0 && adjM[j, i] > 0) {  # Mutualism
            A[i, j] <- A[j, i] <- TRUE
            adjM[i, j] <- runif(1, 0, 2 * adjM[i, j])
            adjM[j, i] <- runif(1, 0, 2 * adjM[j, i])
            Mut[j] <- Mut[i] <- TRUE
          } else if (adjM[i, j] > 0 && adjM[j, i] == 0) {  # Commensalism
            A[i, j] <- TRUE
            adjM[i, j] <- runif(1, 0, 2 * adjM[i, j])
            Bas[j] <- FALSE
          } else if (adjM[i, j] < 0 && adjM[j, i] == 0) {  # Amensalism
            A[i, j] <- TRUE
            adjM[i, j] <- -runif(1, 0, -2 * adjM[i, j])
          }
        }
      }
    }
  }
  
  for (i in 1:rho) {
    if (Mut[i]) {
      adjM[i, i] <- -runif(1, 0, 2 * selfLimAvg[1])
      r[i] <- ifelse(rbinom(1, 1, 0.5) == 1, runif(1, 0, 1), -runif(1, 0, 1))  # Random obligate mutualism
    } else if (Bas[i]) {
      adjM[i, i] <- -runif(1, 0, 2 * selfLimAvg[2])
      r[i] <- runif(1, 0, 1)
    } else {
      adjM[i, i] <- -runif(1, 0, 2 * selfLimAvg[3])
      r[i] <- -runif(1, 0, 1)
    }
    m[i] <- runif(1, 0, 2 * migrAvg)
  }
  
  return(list(interM = adjM, r = r, m = m))
}