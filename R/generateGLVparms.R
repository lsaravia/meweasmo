#' R function to generate random parameters for the Lotka-Volterra model from an adjacency matrix. 
#' 
#' There is an equivalent Rcpp function which presumably has memory issues.
#' 
#' @param adjM       Numeric matrix, metacommunity interaction matrix, the values and the signs are
#'                   the parameters of the LV model then a predator prey interaction will be adjM[i,j]=-1
#'                   adjM[j,i]=1 where j is the predator and i the prey. Competitive interaction is adjM[i,j]=adjM[j,i]=-1
#' @param ef         Efficiency of predator prey interactions
#' @param predIntAvg Average value for the interaction intensity of predator-preys, competition, mutalistic.
#' @param selfLimAvg Numeric vector, and  also set the average value for diagonal entries of the interaction matrix
#'                   that represent self-limitation, the elements of the vector represent 1=mutualistic, 2=Basal, 3=predator species.
#' @param migrAvg    Average value to generate a uniform random migration with range [0,migrAvg*2]
#' @param preserveInt if 0 the values are random defined by `rndType` with `predIntAvg` as a mean, if 1 the values of the interactions adjM[i,i] when i!=j are preserved.
#'                    if 2 the values are random uniform with mean given by adjM[i,i] when i!=j.
#' @param predIntSd  Standard deviation of the distribution to generate the interaction values for `preserveInt = 0` and `rndType = 1` gamma distribution.
#'                   For uniform distribution the sd is fixed at `predIntAvg*2 / sqrt(12)`
#' @param rndType    0=uniform, 1=gamma NOTE: this is implemented only for `preserveInt = 0`
#'
#' @return           A list with the interaction matrix interM, the intrinsic growth rates r, and migration values m
#' 
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
    rate <- predIntSd^2 / predIntAvg  # This is the scale parameter
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