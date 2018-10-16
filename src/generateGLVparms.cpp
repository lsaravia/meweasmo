#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

#define DBG(MSG,X) Rcpp::Rcout << MSG << X << std::endl;
#define DBG4(MSG,X,XX,Y,Z) Rcpp::Rcout << MSG << X << " " << XX << " " << Y << " " << Z << std::endl;
//#define DBG(MSG,X) 
//#define DBG4(MSG,X,XX,Y,Z) 


//’ Generate random parameters for Lotka-Volterra generalized interaction matrix from Lotka-Volterra adjacency matrix  
//’
//’ @param adjM       NumericMatrix, metacommunity interaction matrix, the values and the signs are 
//’                   the parameters of the LV model then a predator prey interaction will be adjM(i,j)=-1 
//’                   adjM(j,i)=1 where j is the predator and i the prey. Competitive interaction is adjM(i,j)=adjM(j,i)=-1
//’ @param ef         Eficiency of predator prey interactions
//’ @param predIntMax Maxium value for the interaction intensity of predator-preys 
//’ @return           A list with the final the number of species by time S, the number of links by time L, the number of basal species
//'                   and the adjacency matrix A. 
// [[Rcpp::export]]
List generateGLVparmsFromAdj(NumericMatrix adjM, double ef, double predIntMax=0.01 ) {
  
  auto rho = adjM.nrow();     // meta web should be a square matrix
  
  LogicalMatrix A(rho,rho);   // Local adyacency matrix
  NumericVector r(rho);     // Growth rates 
  NumericVector m(rho);     // Migration rates
  LogicalVector Bas(rho);    // Basal species

  // Set to all species basal
  std::fill(Bas.begin(), Bas.end(), true);
  
  DBG("Basales: \n",Bas)
  DBG("A: \n",A)
  DBG("adjM: \n", adjM)    
  // Fill adyacency matrix A and number of species SL
  //
  A.fill_diag(true);
  
  for(auto i=0; i<rho; i++ ) {
    for( auto j=0; j<rho; j++) {
      if(adjM(i,j)!=0 && !A(i,j)){
        
        if(adjM(i,j)<0 && adjM(j,i)>0) {  // j is the predator i is the prey
          
          A(i,j)=A(j,i)=true;
          adjM(i,j)=-runif(1,0,predIntMax)[0];         // Random uniform number
          adjM(j,i)= ef*-adjM(i,j);
          Bas[j]=false;                                 // Predators are not basal
          DBG("Predator-Prey ", adjM(i,j))    
            
        } else if(adjM(i,j)<0 && adjM(j,i)<0) { // Competition

          A(i,j)=A(j,i)=true;
          adjM(i,j)=-runif(1,0,predIntMax)[0];
          adjM(j,i)=-runif(1,0,predIntMax)[0];
          DBG("Competition ", adjM(i,j))    
            
        } else if( adjM(i,j)>0 && adjM(j,i)>0 ){        // Mutualism
          
          A(i,j)=A(j,i)=true;
          adjM(i,j)=runif(1,0,predIntMax)[0];
          adjM(j,i)=runif(1,0,predIntMax)[0];
          Bas[j]=false;                                 // Predators are not basal
          DBG("Mutualism ", adjM(i,j))    
            
        } else if( adjM(i,j)>0 && adjM(i,j)==0 ){    // Comensalism
          
          A(i,j)=true;
          adjM(i,j)=runif(1,0,predIntMax)[0];
          Bas[j]=false;
          DBG("Comensalism ", adjM(i,j))    
            
        } else if( adjM(i,j)<0 && adjM(i,j)==0 ){    // Amensalism
          
          A(i,j)=true;
          adjM(i,j)=-runif(1,0,predIntMax)[0];
          DBG("Amensalism ", adjM(i,j))    
            
        }
        DBG("adjM ",adjM)
        DBG("A ", A)
      } 
    }
  }
    
  for(auto i=0;i<rho; i++ ) {
      if(Bas[i]) {
        adjM(i,i) = -runif(1,0,predIntMax*0.1)[0];
        r[i] = runif(1,0,1)[0];
        m[i] = runif(1,0,1)[0];
        
      }
      else {
        adjM(i,i) = -runif(1,0,predIntMax*0.01)[0];
        r[i] = -runif(1,0,1)[0];
        m[i] = runif(1,0,1)[0];
      }
    }
  DBG("Basales: \n",Bas)
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

A1 <- generateGLVparmsFromAdj(A,0.1)


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

A <- matrix(c(0,   1,     0,    0,
              1,   0,     1,    0,
              0,   1,     0,    1,
              0,   0,     1,    0),nrow = 4,byrow=TRUE)


A1 <- generateGLVparms(A,0.1)
#A1$m <- c(0,0,0)
ini <- c(10,10,10,10)
#NumericMatrix metaW,NumericVector m, NumericVector r, NumericVector ini, int time, double tau=0.01)
A2 <- metaWebNetAssemblyGLV(A1$interM,A1$m,A1$r,ini,200,0.1)
A2$STime

A1 <- generateGLVparms(A,0.1)
#A1$m <- c(0,0,0)
ini <- c(10,10,10,10)
#NumericMatrix metaW,NumericVector m, NumericVector r, NumericVector ini, int time, double tau=0.01)
A2 <- metaWebNetAssemblyGLV(A1$interM,A1$m,A1$r,ini,200,0.1)
A2$STime
*/