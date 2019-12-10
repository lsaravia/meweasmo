#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

//#define DBG(MSG,X) Rcpp::Rcout << MSG << X << std::endl;
//#define DBG4(MSG,X,XX,Y,Z) Rcpp::Rcout << MSG << X << " " << XX << " " << Y << " " << Z << std::endl;
#define DBG(MSG,X) 
#define DBG4(MSG,X,XX,Y,Z) 


//' Generate random parameters for Lotka-Volterra generalized interaction matrix from Lotka-Volterra adjacency matrix  
//'
//' @param adjM       NumericMatrix, metacommunity interaction matrix, the values and the signs are 
//'                   the parameters of the LV model then a predator prey interaction will be adjM(i,j)=-1 
//'                   adjM(j,i)=1 where j is the predator and i the prey. Competitive interaction is adjM(i,j)=adjM(j,i)=-1
//' @param ef         Eficiency of predator prey interactions
//' @param predIntAvg Average value for the interaction intensity of predator-preys, competition, mutalistic. 
//' @param selfLimAvg Numeric vector, and  also set the average value for diagonal entries of the interaction matrix
//'                   that represent self-limitation, the elements of the vector represent 1=mutualistic, 2=Basal, 3=predator species.
//' @param migrAvg    Average value to generate at uniform random migration from [0,migrAvg*2] 
//' @param preserveInt if 0 the values are random uniform with prdIntAvg as a mean, if 1 the values of the interactions adjM[i,i] when i!=j are preserved.
//'                    if 2 the values are random uniform with mean given by adjM[i,i] when i!=j.
//' @return           A list with the interaction matrix interM, the intrinsic growth rates r, and migration values m
//' @export
// [[Rcpp::export]]
List generateGLVparmsFromAdj(NumericMatrix adjM, double ef, double predIntAvg=0.01, NumericVector selfLimAvg = NumericVector::create(0.01, 0.01, 0.01),
                             double migrAvg=0.0, int preserveInt=0) {
  
  auto rho = adjM.nrow();     // meta web should be a square matrix
  
  LogicalMatrix A(rho,rho);   // Local adyacency matrix
  NumericVector r(rho);     // Growth rates 
  NumericVector m(rho);     // Migration rates
  LogicalVector Bas(rho);    // Basal species
  LogicalVector Mut(rho);    // Mutualistic species
  
  // Set to all species basal
  std::fill(Bas.begin(), Bas.end(), true);
  std::fill(Mut.begin(), Mut.end(), false);
  
  DBG("Basales: \n",Bas)
  DBG("A: \n",A)
  DBG("adjM: \n", adjM)    
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
            adjM(i,j)=-runif(1,0,predIntAvg*2)[0];         // Random uniform number
            adjM(j,i)= ef*-adjM(i,j);
            Bas[j]=false;                                 // Predators are not basal
            DBG("Predator-Prey ", adjM(i,j))    
              
          } else if(adjM(i,j)<0 && adjM(j,i)<0) { // Competition
  
            A(i,j)=A(j,i)=true;
            adjM(i,j)=-runif(1,0,predIntAvg*2)[0];
            adjM(j,i)=-runif(1,0,predIntAvg*2)[0];
            DBG("Competition ", adjM(i,j))    
              
          } else if( adjM(i,j)>0 && adjM(j,i)>0 ){        // Mutualism
            
            A(i,j)=A(j,i)=true;
            adjM(i,j)=runif(1,0,predIntAvg*2)[0];
            adjM(j,i)=runif(1,0,predIntAvg*2)[0];
            DBG("Mutualism ", adjM(i,j))   
            Mut[j]=Mut[i]=true;
              
          } else if( adjM(i,j)>0 && adjM(j,i)==0 ){    // Comensalism
            
            A(i,j)=true;
            adjM(i,j)=runif(1,0,predIntAvg*2)[0];
            Bas[j]=false;
            DBG("Comensalism ", adjM(i,j))    
              
          } else if( adjM(i,j)<0 && adjM(j,i)==0 ){    // Amensalism
            
            A(i,j)=true;
            adjM(i,j)=-runif(1,0,predIntAvg*2)[0];
            DBG("Amensalism ", adjM(i,j))    
              
          }
          DBG("adjM ",adjM)
          DBG("A ", A)
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
  DBG("Basales: \n",Bas)
  DBG("Mutualistas: \n",Mut)
  DBG("A: \n",A)
  DBG("adjM: \n", adjM)    
    
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

#A1$m <- c(0,0,0)
ini <- c(100,100,100,100)
#NumericMatrix metaW,NumericVector m, NumericVector r, NumericVector ini, int time, double tau=0.01)
A2 <- metaWebNetAssemblyGLV(A1$interM,A1$m,A1$r,ini,200,0.1)
A2$STime


A1 <- generateGLVparmsFromAdj(A,0.1,selfLimAvg=c(0.01,0.01,0.001),preserveInt = 0, migrAvg = 0.2)
ini <- c(0,0,0,0)
A2 <- metaWebNetAssemblyGLV(A1$interM,A1$m,A1$r,ini,200,0.1)
A2$STime


*/