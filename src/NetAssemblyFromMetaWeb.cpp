#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

//#define DBG(MSG,X) Rcpp::Rcout << MSG << X << std::endl;
#define DBG(MSG,X) 




//’ Simulation of an Assembly process from a Meta-web assuming the interactions conserve in the local web 
//’
//’ @param metaW  metacommunity adyacency matrix 
//’ @param m      migration rate (probability) from the meta-web   
//’ @param q      probability of interaction 
//’ @param a      mortality parameter for C   
//’ @param time   Number of time steps of simulation
//’ @return       A list with the final the number of species by time S, the number of links by time L, the number of basal species
//'               and the adjacency matrix A. 
// [[Rcpp::export]]
List metaWebNetAssembly(LogicalMatrix metaW,double m, double q, double a, int time) {
  
  IntegerVector SL(time);      // Species with links vector, the total number of species could be bigger. 

  IntegerVector LL(time);
  auto rho = metaW.nrow();     // meta web should be a square matrix
  
  LogicalMatrix A(rho,rho);   // Local adyacency matrix
  LogicalVector Spc(rho);     // Vector of present species
  LogicalVector Bas(rho);
  unsigned long L=0;
  
  // detect Basal species
  auto r=colSums(metaW);
  for(auto s=0; s<rho; s++){
    if(r[s]==0)
      Bas[s]=true;
  }
  DBG("Basales : ",Bas)
    
  double e = a;  // No extinctions dependent on Connectivity
  
  for(auto t=1; t<time; t++){
    
    DBG("\n============ t: ",t)
    
    // Add Species m*(rho-S)
    //
      
    for(auto s=0; s<rho; s++){
      auto r=runif(1)[0];
        if(r<m)
          if(!Spc[s])
            Spc[s]=true;
      }
    // Update number of species
    // SL[t]=sum(Spc);
    DBG("Spc: ", Spc)  
    DBG("SLt:  ",SL[t])
      
    // Add interactions that exists in the metaweb then j<--i j=predator i=prey
    //  
    for(auto i=0; i<rho; i++ ) {
      for( auto j=0; j<rho; j++) {
        if( Spc[i] && Spc[j] && metaW(i,j)){
            auto r=runif(1)[0];
            if(r<q) {
              A(i,j)=true; // i-->j
            }
          }
        }
      }
      DBG("Interactions A: \n",A)
      L = sum(A);
      DBG("SL ", SL)
        
      // Extinctions  a L / S^2
      //
      //double e=a*L/(S[t]*S[t]);
      //e = rpois(1, e)[0];
      for(auto i=0; i<rho; i++ ) {
        if(Spc[i]){
          auto r=runif(1)[0];
          if(r<e) {
            DBG("SE EXTINGUE i ",i)
            // if species i extinct all column of interactions is deleted
            for( auto j=0; j<rho; j++) {   
              A(j,i)=false; // j-->i
            }
            // if species i extinct all row of interactions is deleted
            for( auto j=0; j<rho; j++) {   
              A(i,j)=false; // j-->i
            }
            Spc[i]=false;
          }
        }
      }
    
    // Extinciones Secundarias
    //
    auto noSecondary=false;
    do{
      noSecondary=false;
      auto r=colSums(A);
      DBG("R:", r)
      for(auto i=0; i<rho; i++ ) {
        if(Spc[i]){
          if(r[i]==0) {
            if(!Bas[i]){                // if it is not basal secondary exctinction
              DBG("EXTINCION SECUNDARIA i ",i)
              for( auto j=0; j<rho; j++) {   // Delete row i (predators), column i (preys) is already 0
                A(i,j)=false; // j-->i
              }
              noSecondary=true;
              Spc[i]=false;
            }
          }
        }
      }
    } while (noSecondary);
    DBG("Extinctions A: \n",A)
    L = sum(A);  // Calculate L after secondary extinctions
    SL[t] = sum((colSums(A)>0)+(rowSums(A)>0)>0); // Species with links
 
    DBG("L :", L)
    DBG("St:",SL[t])
    DBG("Spc:",Spc)
    DBG("Basal:",Bas)
    LL[t]=L;  
    
  }
  unsigned long BB = sum(Bas*Spc);  // Basal species can become exctinct
  
  return List::create(Named("S") = SL,
                      Named("L") = LL,
                      Named("B") = BB,
                      Named("A") = A);

 
}


/*** R
set.seed(123)

# A <- matrix(c(c(0,1,1,1,0),c(0,0,1,0,1),c(0,0,0,0,0),c(0,0,1,1,0),c(0,0,0,1,0)),nrow = 5,byrow=TRUE)

# If q=0 no interactions only basal species survive because secondary extinctions kill all the ones without preys
A0 <- metaWebNetAssembly(A,0.2,0,0,20)
A0

# Cascade Assembly
#
# If m=0 no species 
A0 <- metaWebNetAssembly(A,0,0.2,0,20)
A0

# If q=1 and a=0 a full interaction matrix A is expected, upper triangle less the basal species TRUE
A0 <- metaWebNetAssembly(A,0.2,1,0,20)
A0



# If q=1 and a=0 a full interaction matrix A is expected, upper triangle less the basal species TRUE
A <- matrix(c(c(0,1,1,1,1),c(0,0,1,1,1),c(0,0,0,1,1),c(0,0,0,0,1),c(0,1,0,0,0)),nrow = 5,byrow=TRUE)
A0 <- metaWebNetAssembly(A,0.5,0.5,0.2,20)
A0

# A1 <- cascadeNetAssembly1(20,0.01,0.25,1,10)
# stopifnot(all(A0$S == A0$S))
# stopifnot(all(A0$L == A1$L))
# stopifnot(all(A0$A == A1$A))
# 
# res <- benchmark(cascadeNetAssembly(900,0.01,0.25,1,10),
#                  cascadeNetAssembly1(900,0.01,0.25,1,10),
#                  replications = 30,
#                  order="relative")
# res[,1:4]

# g <- graph_from_adjacency_matrix( AA$A, mode  = "directed")
# dg <- components(g)
# g <- induced_subgraph(g, which(dg$membership == which.max(dg$csize)))
# degree(g,mode="in")
# vcount(g)
# TrophInd(get.adjacency(g,sparse=F))

*/
