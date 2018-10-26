#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

//#define DBG(MSG,X) Rcpp::Rcout << MSG << X << std::endl;
//#define DBG4(MSG,X,XX,Y,Z) Rcpp::Rcout << MSG << X << " " << XX << " " << Y << " " << Z << std::endl;
#define DBG(MSG,X) 
#define DBG4(MSG,X,XX,Y,Z) 



//’ Simulation of an Assembly process from a Meta-web assuming stochastic Generalized Lotka-Volterra
//’ Local dynamics  
//’
//’ @param metaW  NumericMatrix, metacommunity interaction matrix, the values and the signs are the parameters of the LV model  
//’ @param m      NumericVector, migration rate (probability) from the meta-web   
//’ @param r      NumericVector, grow rate of species, could be negative
//’ @param ini    NumericVector, initial condition, or 0  
//’ @param time   Number of time steps of simulation
//’ @param tau    Number of integration steps for Tau leap method 
//’ @return       A list with the final the number of species by time S, the number of links by time L, the number of basal species
//'               and the adjacency matrix A. 
// [[Rcpp::export]]
List metaWebNetAssemblyGLV(NumericMatrix metaW,NumericVector m, NumericVector r, NumericVector ini, int time, double tau=0.01) {
  
  IntegerVector SL(time);      // Species with links vector, the total number of species could be bigger. 

  IntegerVector LL(time);
  auto rho = metaW.nrow();     // meta web should be a square matrix
  
  LogicalMatrix A(rho,rho);   // Local adyacency matrix
  IntegerVector Spc(rho);     // Vector of species populations
  IntegerMatrix ST(rho,time); // Vector of species populations in time

  // Initial values  
  Spc = ini;
  
  // Fill directed adyacency matrix A and number of species SL
  //
  for(auto i=0; i<rho; i++ ) {
    if(Spc[i]>0){
      SL[0]+=1;                 
      ST(i,0)=Spc[i];
    
      for( auto j=0; j<rho; j++) {
        if(Spc[j]>0 && metaW(i,j)!=0 && i!=j){
          if(metaW(i,j)<0 && metaW(j,i)>0){ // predator prey
            A(i,j)=true;
            A(j,i)=false;
          } else if(metaW(i,j)>0 && metaW(j,i)<0) {
            A(i,j)=false;
            A(j,i)=true;
          } else
            A(i,j)=true;
        }
      }
    }
  }
  LL[0] = sum(A);  
  
  
  DBG("Initial values: ",Spc)
  DBG("SL[t]:  ",SL[0])
  DBG("LL[t]:  ",LL[0])
    
  // Calculate the number of inner steps
  auto tauSteps = 1.0/tau;
  
  for(auto t=1; t<time; ++t){
    
    DBG("\n============ t: ",t)
    
    for(auto n=0; n<tauSteps; ++n){

      // Lotka-Volterra dynamics Tau-leap method 
      //  
      for(auto i=0; i<rho; i++ ) {

        auto tauRate =0.0;
        
        for( auto j=0; j<rho; j++) {
            tauRate +=  Spc[j]*metaW(i,j);
            }
        
        DBG4("Spc[i] r tauRate m ",Spc[i],r[i], tauRate,m[i])
        tauRate = ((Spc[i]*(r[i]+tauRate))+m[i])*tau;
        DBG("tauRate ",tauRate)
        if(tauRate>0)
          Spc[i] += rpois(1, tauRate)[0];
        else
          Spc[i] -= rpois(1, -tauRate)[0];
        
        if(Spc[i]<0) Spc[i]=0;
        DBG("Populations: ",Spc)
      }
    }
    
    // End of tauStep cycle

    // Fill adyacency matrix A to determine directed Links and number of species SL 
    //
    std::fill(A.begin(),A.end(),false);
    for(auto i=0; i<rho; i++ ) {
      if(Spc[i]>0){
        SL[t]+=1;    
        ST(i,t)=Spc[i];
        
        for( auto j=0; j<rho; j++) {
          if(Spc[j]>0 && metaW(i,j)!=0 && i!=j){
            if(metaW(i,j)<0 && metaW(j,i)>0){ // predator prey
              A(i,j)=true;
              A(j,i)=false;
            } else if(metaW(i,j)>0 && metaW(j,i)<0) {
              A(i,j)=false;
              A(j,i)=true;
            } else
              A(i,j)=true;
          }
        }
      } 
    }
    
    // number of directed links LL
    //
    LL[t] = sum(A);  

    DBG("Links :", LL[t])
    DBG("Species:",SL[t])
    DBG("Spc:",Spc)

  }

  return List::create(Named("S") = SL,
                      Named("L") = LL,
                      Named("STime") = ST,
                      Named("A") = A);

 
}


/*** R
set.seed(123)

# Predator Prey dynamics
#
A <- matrix(c(-0.001, -0.08,
               0.01, 0   ),nrow = 2,byrow=TRUE)

r <- c(1.5,-0.25)
m <- c(.1,.2)
ini <- c(0,0)
#NumericMatrix metaW,NumericVector m, NumericVector r, NumericVector ini, int time, double tau=0.01)
A1 <- metaWebNetAssemblyGLV(A,m,r,ini,20,0.1)


# 3 Species Chain dynamics
#
A <- matrix(c(-0.001, -0.01, 0,
                0.001,  0   , -0.01,
                   0,  0.002, 0),nrow = 3,byrow=TRUE)

r <- c(0.1,-0.001,-0.03)
m <- c(.1,.1,.1)
ini <- c(10,10,10)
#NumericMatrix metaW,NumericVector m, NumericVector r, NumericVector ini, int time, double tau=0.01)
A1 <- metaWebNetAssemblyGLV(A,m,r,ini,100,0.1)
A1
df <- as.data.frame(t(A1$STime))
df$Time <- 1:100
require(tidyr)
df <- gather(df,key="Species",value="N", V1:V3)
require(ggplot2)
ggplot(df, aes(Time,N,colour=Species)) + geom_line()


*/
