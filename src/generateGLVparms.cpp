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
//’                   the parameters of the LV model  
//’ @param ef         Eficiency of predator prey interactions
//’ @param predIntMax Maxium value for the interaction intensity of predator-preys 
//’ @return           A list with the final the number of species by time S, the number of links by time L, the number of basal species
//'                   and the adjacency matrix A. 
// [[Rcpp::export]]
List generateGLVparms(NumericMatrix adjM, double ef, double predIntMax=0.01 ) {
  
  auto rho = adjM.nrow();     // meta web should be a square matrix
  
  LogicalMatrix A(rho,rho);   // Local adyacency matrix
  NumericVector r(rho);     // Growth rates 
  NumericVector m(rho);     // Migration rates
  LogicalMatrix Bas(rho);    // Basal species

  Bas = true;               // Set to all species basal
  
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
          
        } else if(adjM(i,j)<0 && adjM(j,i)<0) { // Competition

          A(i,j)=A(j,i)=true;
          adjM(i,j)=-runif(1,0,predIntMax)[0];
          adjM(j,i)=-runif(1,0,predIntMax)[0];
          
        } else if( adjM(i,j)>0 && adjM(i,j)>0 ){        // Mutualism
          
          A(i,j)=A(j,i)=true;
          adjM(i,j)=runif(1,0,predIntMax)[0];
          adjM(j,i)=runif(1,0,predIntMax)[0];
          
        } else if( adjM(i,j)>0 && adjM(i,j)==0 ){    // Comensalism
          
          A(i,j)=true;
          adjM(i,j)=runif(1,0,predIntMax)[0];

        } else if( adjM(i,j)<0 && adjM(i,j)==0 ){    // Amensalism
          
          A(i,j)=true;
          adjM(i,j)=-runif(1,0,predIntMax)[0];
          
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
        m[i] = runif(1,0,predIntMax*0.01)[0];
        
      }
      else {
        adjM(i,i) = -runif(1,0,predIntMax*0.01)[0];
        r[i] = -runif(1,0,1)[0];
        m[i] = runif(1,0,predIntMax*0.01)[0];
      }
    }
  DBG("Basales: \n",Bas)
  DBG("A: \n",A)
  DBG("adjM: \n", adjM)    
    
  return List::create(Named("interM") = adjM,
                      Named("r") = r,
                      Named("m") = m);
  
  
}