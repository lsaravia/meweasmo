#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

//#define DBG(MSG,X) Rcpp::Rcout << MSG << X << std::endl;
//#define DBG4(MSG,X,XX,Y,Z) Rcpp::Rcout << MSG << X << " " << XX << " " << Y << " " << Z << std::endl;
#define DBG(MSG,X) 
#define DBG4(MSG,X,XX,Y,Z) 


//’ Generate random Lotka-Volterra adjacency matrix with fixed proportion of interactions   
//’
//’ @param numSp      Integer, number of species 
//’ @param C          Double, Probability of interaction = Connectivity = number of interacions/(numSp^2)
//’ @param propInt    Numeric vector, Proportion of interactions, there are 5 types of interactions the elements of the vector
//’                   represent 1= mutualism (+,+), 2=Predator-prey (+,-), 3= competition (-,-), 4=Comensalims (+,0), 
//’                   5=Amensalism (-,0). 
//’ @return           A NumericMatrix with the GLV type of adjacency matrix  
// [[Rcpp::export]]
IntegerMatrix generateRandomGLVadjMat(int numSp, double C, NumericVector propInt) {
  
  auto rho = propInt.length();     // meta web should be a square matrix
  if(rho!=5) 
    throw std::range_error("propInt length different from 5");
  
  IntegerMatrix A(numSp,numSp);      // Local adyacency matrix
  IntegerMatrix adjM(numSp,numSp);   // Local adyacency matrix
  
  for(auto i=0; i<numSp; i++ ) {
      for( auto j=0; j<numSp; j++) { 
        if( i!=j){ 
         A(i,j) = rbinom(1,1,C)[0];
        }
      }
  }
  NumericVector accum(5);
  std::partial_sum(propInt.begin(), propInt.end(), accum.begin());
  //accum = cumsum(propInt);
  
  DBG("adjM: \n", adjM)    
  DBG("accum: \n", accum)    
    
  // Fill adyacency matrix A and number of species SL
  //
  for(auto i=0; i<numSp; i++ ) {
    for( auto j=0; j<numSp; j++) {
      if(adjM(i,j)==0 && adjM(j,i)==0 && A(i,j)){
        
        DBG4("A(i,j) ", A(i,j),"i,j",i,j)
        
        auto s = runif(1,0,1)[0];
        auto p = 0;

        DBG("s: ",s)
          
        while(s>accum[p]){
          p++;
        }
        
        DBG("p: ",p)
          
        switch(p){
        case 0:
          // Mutualistic
          adjM(i,j)=adjM(j,i)=1;
          break;
        
        case 1:
          // Predator-Prey j is the predator
          adjM(i,j)=-1;
          adjM(j,i)=1;
          break;
          
        case 2:
          // Competition 
          adjM(i,j)=adjM(j,i)=-1;
          break;
            
        case 3:
          // Comensalism j beneficts i
          adjM(i,j)=1;
          adjM(j,i)=0;
          break;
        
        case 4:
          // Amensalism  j harms i
          adjM(i,j)=-1;
          adjM(j,i)=0;
          break;
       }
      }
    }
  }
    

  DBG("A: \n",A)
  DBG("adjM: \n",adjM)
        
  return adjM;
  
  
}




/*** R
set.seed(123)

# Predator Prey dynamics
#

# Generate a network with all predator-prey (antagonistic) and some basal species 
A <- generateRandomGLVadjMat(4,0.5,c(0.0,.8,0,0,0)) 
A
A1 <- generateGLVparmsFromAdj(A,0.1)
A1

# Generate a network with all mutualistic  
A <- generateRandomGLVadjMat(4,0.5,c(.8,0,0,0,0)) 
A
A1 <- generateGLVparmsFromAdj(A,0.1)
A1

# Generate a network with all competitive  
A <- generateRandomGLVadjMat(4,0.5,c(0,0,1,0,0)) 
A
A1 <- generateGLVparmsFromAdj(A,0.1)
A1

# Generate a network with all comensalism  
A <- generateRandomGLVadjMat(4,0.5,c(0,0,0,1,0)) 
A
A1 <- generateGLVparmsFromAdj(A,0.1)
A1

# Generate a network with all amensalism  
A <- generateRandomGLVadjMat(4,0.5,c(0,0,0,0,1)) 
A
A1 <- generateGLVparmsFromAdj(A,0.1)
A1


# Generate a network with all types  
A <- generateRandomGLVadjMat(10,0.5,c(0.2,0.2,0.2,0.2,0.2)) 
A
A1 <- generateGLVparmsFromAdj(A,0.1,0.01)
A1

#A1$m <- c(0,0,0,0)
ini <- rep(0,times=nrow(A1$interM))

#NumericMatrix metaW,NumericVector m, NumericVector r, NumericVector ini, int time, double tau=0.01)
A2 <- metaWebNetAssemblyGLV(A1$interM,A1$m,A1$r,ini,200,0.1)
A2$STime
A2$S
*/