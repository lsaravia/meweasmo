#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

//#define DBG(MSG,X) Rcpp::Rcout << MSG << X << std::endl;
//#define DBG4(MSG,X,XX,Y,Z) Rcpp::Rcout << MSG << X << " " << XX << " " << Y << " " << Z << std::endl;
#define DBG(MSG,X) 
#define DBG4(MSG,X,XX,Y,Z) 


//' Calculate the proportion of interacions types in a Lotka-Volterra adjacency matrix
//'
//' @param adjM       Numeric Matrix Lotka-Volterra adjcency matrix  
//' @param sp         Numeric vector with the present species, if species=0 the interaction is not counted
//' 
//' @return           Numeric vector with 5 elements Matrix 1=Mutualistic (+,+), 2=Predator-prey (-,+), 3=Competitive (-,-)
//'                   4=Comensalism (+,0), 5=Amensalism (-,0)
//' @export
// [[Rcpp::export]]
NumericVector calcPropInteractionsGLVadjMat(NumericMatrix adjM, NumericVector spc) {

  auto numSp = adjM.nrow();     // meta web should be a square matrix
  
  NumericMatrix A(clone(adjM));      // Local adyacency matrix
  
  NumericVector propInt(5);
  DBG("adjM: \n", adjM)    
  DBG("propInt: \n", propInt)    
    
  // Fill adyacency matrix A and number of species SL
  //
  for(auto i=0; i<numSp; i++ ) {
    for( auto j=0; j<numSp; j++) {
      if(i!=j && spc[i]>0 && spc[j]>0){
        
        DBG4("adjM(i,j) ", adjM(i,j),"i,j",i,j)
        if(A(i,j)>0 && A(j,i)>0 ){
          // Mutualistic
          propInt[0]+=2;
          A(i,j)=A(j,i)=0;
        }
        else if(A(i,j)<0 && A(j,i)>0 ){
          // Predator-Prey j is the predator
          propInt[1]++;
          A(i,j)=A(j,i)=0;
        }
        else if(A(i,j)>0 && A(j,i)<0 ){
          // Predator-Prey i is the predator
          propInt[1]++;
          A(i,j)=A(j,i)=0;
        }         
        else if(A(i,j)<0 && A(j,i)<0 ){
          // Competition 
          propInt[2]+=2;
          A(i,j)=A(j,i)=0;
        }
        else if(A(i,j)>0 && A(j,i)==0 ){
          // Comensalism j beneficts i
          propInt[3]++;
          A(i,j)=0;
        }
        else if(A(i,j)<0 && A(j,i)==0 ){
          // Amensalism  j harms i
          propInt[4]++;
          A(i,j)=0;
        }        
       }
      }
    }
  auto tot =  sum(propInt);
  DBG("propInt: \n", propInt);
  DBG("tot: ", tot);
  if(tot>0)
    propInt = propInt/tot;
  return( propInt);
  
}
    




/*** R
set.seed(123)


# Generate a network with all types  
A <- generateRandomGLVadjMat(5,0.5,c(0.2,0.2,0.2,0.2,0.2)) 
A[3,1] <- 1
A1 <- generateGLVparmsFromAdj(A,0.1,0.01)
A1
ini <- rep(1,times=nrow(A1$interM))

calcPropInteractionsGLVadjMat(A1$interM,ini)



#NumericMatrix metaW,NumericVector m, NumericVector r, NumericVector ini, int time, double tau=0.01)
A2 <- metaWebNetAssemblyGLV(A1$interM,A1$m,A1$r,ini,200,0.1)
A2$STime
A2$S
A2$L
A2$L/(A2$S*A2$S)

calcPropInteractionsGLVadjMat(A1$interM,A2$STime[,200])
A2$A


# 3 Species Chain dynamics
#
A <- matrix(c(-0.001, -0.01, 0,
              0.001,  0   , -0.01,
              0,  0.002, 0),nrow = 3,byrow=TRUE)
A1 <- generateGLVparmsFromAdj(A,0.1,0.01)
A1
ini <- rep(1,times=nrow(A1$interM))

calcPropInteractionsGLVadjMat(A1$interM,ini)

*/