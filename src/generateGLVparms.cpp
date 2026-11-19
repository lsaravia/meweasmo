#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

//#define DBG(MSG,X) Rcpp::Rcout << MSG << X << std::endl;
//#define DBG4(MSG,X,XX,Y,Z) Rcpp::Rcout << MSG << X << " " << XX << " " << Y << " " << Z << std::endl;
#define DBG(MSG,X) 
#define DBG4(MSG,X,XX,Y,Z) 


//' Generate random parameters of a generalized Lotka-Volterra model from an adjacency matrix
//'
//' Builds the parameters needed by [metaWebNetAssemblyGLV()] (interaction matrix, intrinsic growth rates 
//' and migration rates) from a matrix that only defines the structure of the interactions, 
//' for example one generated with [generateRandomGLVadjMat()]. The input matrix is not modified.
//' 
//' The type of each interaction is read from the signs of `adjM` (`adjM[i, j]` is the effect of 
//' species `j` on species `i`):
//' * `adjM[i, j] < 0` and `adjM[j, i] > 0`: predator-prey, `j` is the predator and `i` the prey. 
//'   The strength of the effect on the prey is random and the effect on the predator is 
//'   `ef` times the effect on the prey, with opposite sign.
//' * both negative: competition; both positive: mutualism. Each direction gets its own random strength.
//' * `adjM[i, j] > 0` and `adjM[j, i] == 0`: commensalism (`i` benefits); 
//'   `adjM[i, j] < 0` and `adjM[j, i] == 0`: amensalism (`i` is harmed).
//' 
//' The strengths of the interactions (off-diagonal values) are drawn according to `preserveInt` and `rndType`.
//' 
//' The remaining parameters depend on the type of each species, which is also deduced from `adjM`: 
//' *mutualistic* if it takes part in any mutualistic interaction, otherwise *basal* if it has no
//' resources (it is never the predator in a predator-prey interaction, nor the species `j` in a commensalism), 
//' otherwise *consumer*. For each type:
//' * Self-limitation (diagonal of the interaction matrix): `-runif(0, 2 * selfLimAvg[k])`, 
//'   where `k` is 1 for mutualistic, 2 for basal and 3 for consumer species. 
//'   The diagonal of `adjM` is ignored.
//' * Intrinsic growth rate `r`: positive, `runif(0, 1)`, for basal species; negative, `-runif(0, 1)`, for 
//'   consumers (they need resources to persist); and of random sign, with the same magnitude, for 
//'   mutualistic species (random obligate mutualism).
//' * Migration `m`: `runif(0, 2 * migrAvg)` for all species.
//' 
//' There is an equivalent implementation in R, [generateGLVparmsFromAdjR()], that gives the same 
//' results for the same random seed.
//'
//' @param adjM       Numeric square matrix with the structure of the interactions. Only the signs and the zeros
//'                   are used when `preserveInt = 0`; with `preserveInt` 1 or 2 the values of the non-zero 
//'                   off-diagonal entries are used as well. 
//' @param ef         Efficiency of predator-prey interactions: effect on the predator / effect on the prey.
//' @param predIntAvg Average strength of the interactions (predator-prey, competition and mutualism) 
//'                   when `preserveInt = 0`.
//' @param selfLimAvg Numeric vector of length 3 with the average self-limitation of 
//'                   1 = mutualistic, 2 = basal and 3 = consumer species (diagonal entries of the interaction matrix).
//' @param migrAvg    Average migration rate, the migration of each species is uniform in \[0, 2 * `migrAvg`\].
//' @param preserveInt How to define the strength of the interactions:
//' * 0: random values defined by `rndType` with mean `predIntAvg`.
//' * 1: the values of `adjM` are preserved (off-diagonal entries only).
//' * 2: random uniform values with mean given by the value in `adjM` for each interaction, 
//'   that is, uniform in \[0, 2 * abs(`adjM[i, j]`)\].
//' @param predIntSd  Standard deviation of the strengths when `preserveInt = 0` and `rndType = 1` (gamma distribution).
//'                   For the uniform distribution the sd is fixed at `predIntAvg * 2 / sqrt(12)`.
//' @param rndType    Distribution of the strengths when `preserveInt = 0`: 0 = uniform, 1 = gamma. 
//'                   Other values of `preserveInt` always use the uniform distribution.
//' 
//' @return A list with the interaction matrix `interM` (with the self-limitation in the diagonal), 
//'         the intrinsic growth rates `r` and the migration rates `m`. They can be passed to 
//'         [metaWebNetAssemblyGLV()] as `metaW`, `r` and `m`.
//' @seealso [generateRandomGLVadjMat()] to generate `adjM`, [generateGLVparmsFromAdjR()] for the R 
//'          implementation, and [metaWebNetAssemblyGLV()].
//' @examples
//' set.seed(1)
//' adj <- generateRandomGLVadjMat(20, 0.2, c(0.1, 0.3, 0.3, 0.15, 0.15))
//' parms <- generateGLVparmsFromAdj(adj, ef = 0.1, predIntAvg = 0.01, migrAvg = 0.01)
//' str(parms)
//' 
//' # Gamma distributed interaction strengths
//' parms <- generateGLVparmsFromAdj(adj, ef = 0.1, predIntAvg = 0.01, predIntSd = 0.001, rndType = 1)
//' @export
// [[Rcpp::export]]
List generateGLVparmsFromAdj(NumericMatrix adjM, double ef, double predIntAvg=0.01, NumericVector selfLimAvg = NumericVector::create(0.01, 0.01, 0.01),
                             double migrAvg=0.0, int preserveInt=0, double predIntSd=0.01, int rndType=0) {
  
  adjM = clone(adjM);         // do not modify the caller's matrix in place
  auto rho = adjM.nrow();     // meta web should be a square matrix

  LogicalMatrix A(rho,rho);   // Local adyacency matrix
  NumericVector r(rho);     // Growth rates
  NumericVector m(rho);     // Migration rates
  LogicalVector Bas(rho);    // Basal species
  LogicalVector Mut(rho);    // Mutualistic species
  
  // Set to all species basal
  std::fill(Bas.begin(), Bas.end(), true);
  std::fill(Mut.begin(), Mut.end(), false);
  
  auto a_predInt = 0.0, b_predInt = 0.0, shape = 0.0, rate = 0.0;
  if(rndType==0) {
    // Uniform distribution
    //
    // Parameters a and b, sd fixed at predIntAvg*2/sqrt(12)
    a_predInt = 0.0;
    b_predInt = predIntAvg * 2;
  } else if(rndType==1) {
    // Gamma distribution
    //
    DBG("Gamma: \n",rndType)
    shape = pow(predIntAvg/predIntSd,2);
    rate  = pow(predIntSd,2)/predIntAvg;   // This  the scale parameter
  }
    
  //DBG("Basales: \n",Bas)
  //DBG("A: \n",A)
  //DBG("adjM: \n", adjM)    
  
  // Fill adyacency matrix A and number of species SL
  //
  A.fill_diag(true);
  if(preserveInt==0)
  {
    for(auto i=0; i<rho; i++ ) {
      for( auto j=0; j<rho; j++) {
        if(adjM(i,j)!=0 && !A(i,j)){
          
          if(adjM(i,j)<0 && adjM(j,i)>0) {  // j is the predator i is the prey
            
            A(i,j)=A(j,i)=true;
            if(rndType==0)
              adjM(i,j)=-runif(1,a_predInt,b_predInt)[0];         // Random uniform number
            else if(rndType==1)
            {
              adjM(i,j)=-rgamma(1,shape,rate)[0];         // Random gamma number
              DBG("\t", adjM(i,j))
            }
            adjM(j,i)= ef*-adjM(i,j);
            Bas[j]=false;                                 // Predators are not basal
            //DBG("Predator-Prey ", adjM(i,j))    
              
          } else if(adjM(i,j)<0 && adjM(j,i)<0) { // Competition
  
            A(i,j)=A(j,i)=true;
            if(rndType==0) {
              adjM(i,j)=-runif(1,a_predInt,b_predInt)[0];         // Random uniform number
              adjM(j,i)=-runif(1,a_predInt,b_predInt)[0];         // Random uniform number
            } else if(rndType==1) {
              adjM(i,j)=-rgamma(1,shape,rate)[0];         // Random gamma number
              adjM(j,i)=-rgamma(1,shape,rate)[0];         // Random gamma number
              DBG("\t", adjM(j,i))
              DBG("\t", adjM(i,j))
            }
            //DBG("Competition ", adjM(i,j))    
              
          } else if( adjM(i,j)>0 && adjM(j,i)>0 ){        // Mutualism
            
            A(i,j)=A(j,i)=true;
            if(rndType==0) {
              adjM(i,j)=runif(1,a_predInt,b_predInt)[0];         // Random uniform number
              adjM(j,i)=runif(1,a_predInt,b_predInt)[0];         // Random uniform number
            } else if(rndType==1) {
              adjM(i,j)=rgamma(1,shape,rate)[0];         // Random gamma number
              adjM(j,i)=rgamma(1,shape,rate)[0];         // Random gamma number
              DBG("\t", adjM(j,i))
              DBG("\t", adjM(i,j))
            }
            //DBG("Mutualism ", adjM(i,j))   
            Mut[j]=Mut[i]=true;
              
          } else if( adjM(i,j)>0 && adjM(j,i)==0 ){    // Comensalism
            
            A(i,j)=true;
            if(rndType==0)
              adjM(i,j)=runif(1,a_predInt,b_predInt)[0];         // Random uniform number
            else if(rndType==1) {
              adjM(i,j)=rgamma(1,shape,rate)[0];         // Random gamma number
              DBG("\t", adjM(i,j))
            }
            Bas[j]=false;
            //DBG("Comensalism ", adjM(i,j))    
              
          } else if( adjM(i,j)<0 && adjM(j,i)==0 ){    // Amensalism
            
            A(i,j)=true;
            if(rndType==0)
              adjM(i,j)=-runif(1,a_predInt,b_predInt)[0];         // Random uniform number
            else if(rndType==1) {
              adjM(i,j)=-rgamma(1,shape,rate)[0];         // Random gamma number
              DBG("\t", adjM(i,j))
            }
            //DBG("Amensalism ", adjM(i,j))    
              
          }
          //DBG("adjM ",adjM)
          //DBG("A ", A)
        } 
      }
    }
  }
  else if( preserveInt ==1)     // Do not touch the interaction values
  {
    for(auto i=0; i<rho; i++ ) {
      for( auto j=0; j<rho; j++) {
        if(adjM(i,j)!=0 && !A(i,j)){
          
          if(adjM(i,j)<0 && adjM(j,i)>0) {  // j is the predator i is the prey
            
            A(i,j)=A(j,i)=true;
            //adjM(j,i)= ef*-adjM(i,j);
            Bas[j]=false;                                 // Predators are not basal
            DBG("Predator-Prey ", adjM(i,j))    
              
          } else if(adjM(i,j)<0 && adjM(j,i)<0) { // Competition
            
            A(i,j)=A(j,i)=true;
            DBG("Competition ", adjM(i,j))    
              
          } else if( adjM(i,j)>0 && adjM(j,i)>0 ){        // Mutualism
            
            A(i,j)=A(j,i)=true;
            DBG("Mutualism ", adjM(i,j))   
              Mut[j]=Mut[i]=true;
            
          } else if( adjM(i,j)>0 && adjM(j,i)==0 ){    // Comensalism
            
            A(i,j)=true;
            Bas[j]=false;
            DBG("Comensalism ", adjM(i,j))    
              
          } else if( adjM(i,j)<0 && adjM(j,i)==0 ){    // Amensalism
            
            A(i,j)=true;
            DBG("Amensalism ", adjM(i,j))    
              
          }
          DBG("adjM ",adjM)
            DBG("A ", A)
        } 
      }
    }
  }
  else if( preserveInt ==2)
  {
    for(auto i=0; i<rho; i++ ) {
      for( auto j=0; j<rho; j++) {
        if(adjM(i,j)!=0 && !A(i,j)){
          
          if(adjM(i,j)<0 && adjM(j,i)>0) {  // j is the predator i is the prey
            
            A(i,j)=A(j,i)=true;
            DBG("Predator-Prey before runif", adjM(i,j))    
            adjM(i,j)=-runif(1,0,-adjM(i,j)*2)[0];         // Random uniform number
            adjM(j,i)= ef*-adjM(i,j);
            Bas[j]=false;                                 // Predators are not basal
            DBG("Predator-Prey ", adjM(i,j))    
              
          } else if(adjM(i,j)<0 && adjM(j,i)<0) { // Competition
  
            DBG("Competition before runif", adjM(i,j))    
            A(i,j)=A(j,i)=true;
            adjM(i,j)=-runif(1,0,-adjM(i,j)*2)[0];
            adjM(j,i)=-runif(1,0,-adjM(j,i)*2)[0];
            DBG("Competition ", adjM(i,j))    
              
          } else if( adjM(i,j)>0 && adjM(j,i)>0 ){        // Mutualism
            DBG("Mutualism before runif", adjM(i,j))    
            A(i,j)=A(j,i)=true;
            adjM(i,j)=runif(1,0,adjM(i,j)*2)[0];
            adjM(j,i)=runif(1,0,adjM(j,i)*2)[0];
            DBG("Mutualism ", adjM(i,j))   
            Mut[j]=Mut[i]=true;
              
          } else if( adjM(i,j)>0 && adjM(j,i)==0 ){    // Comensalism
            DBG("Comensalism before runif", adjM(i,j))    
            A(i,j)=true;
            adjM(i,j)=runif(1,0,adjM(i,j)*2)[0];
            Bas[j]=false;
            DBG("Comensalism ", adjM(i,j))    
              
          } else if( adjM(i,j)<0 && adjM(j,i)==0 ){    // Amensalism
            
            DBG("Amensalism before runif", adjM(i,j))    
            A(i,j)=true;
            adjM(i,j)=-runif(1,0,-adjM(i,j)*2)[0];
            DBG("Amensalism ", adjM(i,j))    
              
          }
          DBG("adjM ",adjM)
          DBG("A ", A)
        } 
      }
    }
  }
    
    
  for(auto i=0;i<rho; i++ ) {
    if (Mut[i]) {
      adjM(i,i) = -runif(1,0,selfLimAvg[0]*2)[0];
      if(rbinom(1,1,0.5)[0]==1)                                 // Random Obligate mutualism
        r[i] = runif(1,0,1)[0];
      else
        r[i] = -runif(1,0,1)[0];
    }
    else if(Bas[i]) {
        adjM(i,i) = -runif(1,0,selfLimAvg[1]*2)[0];
        r[i] = runif(1,0,1)[0];

    }  
    else {
        adjM(i,i) = -runif(1,0,selfLimAvg[2]*2)[0];
        r[i] = -runif(1,0,1)[0];
      }
    m[i] = runif(1,0,migrAvg*2)[0];
    
    }
  //DBG("Basales: \n",Bas)
  //DBG("Mutualistas: \n",Mut)
  //DBG("A: \n",A)
  //DBG("adjM: \n", adjM)    
    
  return List::create(Named("interM") = adjM,
                      Named("r") = r,
                      Named("m") = m);
  
  
}




/*** R
set.seed(123)

# Predator Prey dynamics
#
A <- matrix(c(-0.001, -0.08,
                0.01, 0   ),nrow = 2,byrow=TRUE)

A1 <- generateGLVparmsFromAdj(A,0.1,selfLimAvg = c(1,0.001,0.01),migrAvg = 0.8)

A <- matrix(c(-0.001, -0.08,
              0.01, 0   ),nrow = 2,byrow=TRUE)

A1 <- generateGLVparmsFromAdj(A,0.1,selfLimAvg = c(1,0.001,0.01),migrAvg = 0.8,preserveInt = 1)

A1

A <- matrix(c(-0.001, -0.08,
              0.01, 0   ),nrow = 2,byrow=TRUE)

A1 <- generateGLVparmsFromAdj(A,0.1,selfLimAvg = c(1,0.001,0.01),migrAvg = 0.8,preserveInt = 2)

A <- matrix(c(-0.001,  -0.01,     0,
               0.001,      0,  -0.01,
                   0,  0.002,     0),nrow = 3,byrow=TRUE)

A <- matrix(c(-0.001,  -0.01,     0,
               0.001,      0,  0.01,
                   0,  -0.002,     0),nrow = 3,byrow=TRUE)

A <- matrix(c(-1,  -1,     0,    0,
               1,   0,     1,    0,
               0,  -1,     0,    1,
               0,   0,     1,    0),nrow = 4,byrow=TRUE)
A1 <- generateGLVparmsFromAdj(A,0.1)

# Mutualism + competition
#
A <- matrix(c(0,   1,     1,    0,
              1,   0,     1,    0,
              1,   1,     0,    -1,
              0,   0,     -1,    0),nrow = 4,byrow=TRUE)

# Test No migration , 0 initial conditions, the results should be always 0
#
A1 <- generateGLVparmsFromAdj(A,0.1,selfLimAvg=c(0.1,0.01,0.001))
#A1$m <- c(0,0,0,0)
ini <- c(0,0,0,0)

#NumericMatrix metaW,NumericVector m, NumericVector r, NumericVector ini, int time, double tau=0.01)
A2 <- metaWebNetAssemblyGLV(A1$interM,A1$m,A1$r,ini,200,0.1)
A2$STime
d <- apply(A2$STime[,101:200],1,mean)>5
A <- A2$A[d,d]
sum(d)
sum(A)

A <- matrix(c(0,   1,     1,    0,
              1,   0,     1,    0,
              1,   1,     0,    -1,
              0,   0,     -1,    0),nrow = 4,byrow=TRUE)

#A1$m <- c(0,0,0)
ini <- c(100,100,100,100)
#NumericMatrix metaW,NumericVector m, NumericVector r, NumericVector ini, int time, double tau=0.01)
A2 <- metaWebNetAssemblyGLV(A1$interM,A1$m,A1$r,ini,200,0.1)
A2$STime


A1 <- generateGLVparmsFromAdj(A,0.1,selfLimAvg=c(0.01,0.01,0.001),preserveInt = 0, migrAvg = 0.2)
ini <- c(0,0,0,0)
A2 <- metaWebNetAssemblyGLV(A1$interM,A1$m,A1$r,ini,200,0.1)
A2$STime


# Generate a mutualistic / Predator-prey network   
A <- generateRandomGLVadjMat(100,0.5,c(0.2,0.2,0.2,0.2,0.2)) 
A <- generateRandomGLVadjMat(100,0.5,c(1,0.2,0.2,0.2,0.2)) 

A
A1 <- generateGLVparmsFromAdj(A,0.1,0.001,selfLimAvg=c(0.01,0.01,0.001),preserveInt = 0, migrAvg = 0.2)
ini <- rep(0.0,100)
A2 <- metaWebNetAssemblyGLV(A1$interM,A1$m,A1$r,ini,200,0.1)
A2$STime[,200]

DG <-A1$interM
diag(DG) <- 0
mean( abs(DG[DG!=0] ))
  */