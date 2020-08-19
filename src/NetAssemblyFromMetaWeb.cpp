#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// #define DBG(MSG,X) Rcpp::Rcout << MSG << X << std::endl;
// #define DBG3(MSG,X,Y,Z) Rcpp::Rcout << MSG << X << "\t" << Y << "\t" << Z << std::endl;
#define DBG(MSG,X)
#define DBG3(MSG,X,Y,Z)




//' Simulation of an Assembly process from a Meta-web assuming the interactions are conserved in the local web 
//'
//' This is a dynamical model of colonization and extinction process, with the restriction 
//' that predators must have at least one prey species to survive in the local community,  described in [1] 
//' With an additional secondary extinctions probability that controls when a predator goes extinct if it has
//' no prey.
//' 
//' @references
//' 
//' 1. Galiana, N., Lurgi, M., Claramunt-López, B., Fortin, M.-J., Leroux, S., Cazelles, K., et al. (2018). 
//' The spatial scaling of species interaction networks. Nat. Ecol. Evol., 2, 782–790
//'
//' @param metaW  metacommunity adyacency matrix 
//' @param m      A numeric vector of species' migration rates (probability) from the meta-web   
//' @param e      A numeric vector of species' extinction probability 
//' @param se      A numeric vector of species' secondary extinction probability 
//' @param time   Number of time steps of simulation
//' @return       A list with the final the number of species by time S, the number of links by time L, 
//'               the time series of species STime and the adjacency matrix A with effective links. 
// [[Rcpp::export]]
List metaWebNetAssembly(LogicalMatrix metaW, NumericVector m, NumericVector e, NumericVector se, int time) {

  if( metaW.ncol()!=metaW.ncol())  
    stop("Matrix metaW have to be square");
  
  IntegerVector SL(time);      // Species with links vector, the total number of species could be bigger. 
  IntegerVector LL(time);
  
  auto rho = metaW.nrow();     // meta web should be a square matrix
  if(m.length()!=rho)
    stop("migration vector m must have the same elements that the dimensions of metaW");  
  if(e.length()!=rho)
    stop("extinction vector e must have the same elements that the dimensions of metaW");  
  
  LogicalMatrix A(rho,rho);   // Local adyacency matrix
  LogicalVector Spc(rho);     // Vector of present species
  LogicalVector Bas(rho);     // Vector of basal species
  IntegerMatrix ST(rho,time); // Vector of species populations in time
  unsigned long L=0;
  
  // detect Basal species
  auto r=colSums(metaW);
  for(auto s=0; s<rho; s++){
    if(r[s]==0)
      Bas[s]=true;
  }
  DBG("Basales : ",Bas)
    
  for(auto t=1; t<time; t++){
    
    DBG("\n============ t: ",t)
    
    // Add Species m*(rho-S)
    //
    auto r=runif(rho);  
    DBG("\nrandom: ",r)
    DBG("\nmigrat: ",m)
      
    for(auto s=0; s<rho; s++){
        if(r[s]<m[s])
          if(!Spc[s]){
            Spc[s]=true;
          }
      }
    // Update number of species
    ST(_,0)=Spc;
    
    // SL[t]=sum(Spc);
    DBG("Spc: ", Spc)  
    DBG("SLt:  ",SL[t])

    // Add interactions that exists in the metaweb then j<--i j=predator i=prey
    //  
    for(auto i=0; i<rho; i++ ) {
      for( auto j=0; j<rho; j++) {
        if( Spc[i] && Spc[j] && metaW(i,j)){
            A(i,j)=true; // i-->j
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
        if(r<e[i]) {
          DBG3("SE EXTINGUE i, rnd, e ",i,r,e[i])
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
            if(!Bas[i]){ 
              // if it is not basal secondary exctinction
              auto rnd = runif(1)[0];
                if(rnd < se[i]) {
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
      }
    } while (noSecondary);
    DBG("Extinctions A: \n",A)
    L = sum(A);  // Calculate L after secondary extinctions
    SL[t] = sum((colSums(A)>0)+(rowSums(A)>0)>0); // Active links
    
    ST(_,t) = Spc;
    
    DBG("L :", L)
    DBG("St:",SL[t])
    DBG("Spc:",Spc)
    DBG("Basal:",Bas)
    LL[t]=L;  
    
  }
  unsigned long BB = sum(Bas*Spc);  // Basal species can become exctinct
  
  return List::create(Named("S") = SL,
                      Named("L") = LL,
                      Named("STime") = ST,
                      Named("A") = A);

 
}


//' Simulation of an Assembly process from a Meta-web assuming the interactions are conserved in the local web 
//'
//' This is a dynamical model of colonization and extinction process, with the restriction 
//' that predators must have at least one prey species to survive in the local community, described in [1]. 
//' With an additional secondary extinctions probability that controls when a predator goes extinct if it has
//' no prey. This is a continuous time version of the model that follows the Gillespie algorithm [2] for simulation 
//' 
//' @references
//' 
//' 1. Galiana, N., Lurgi, M., Claramunt-López, B., Fortin, M.-J., Leroux, S., Cazelles, K., et al. (2018). 
//' The spatial scaling of species interaction networks. Nat. Ecol. Evol., 2, 782–790
//' 2. Gillespie, D. T. (1976). A general method for numerically simulating the stochastic time evolution of coupled chemical reactions. 
//' Journal of Computational Physics, 22(4), 403–434. doi: 10.1016/0021-9991(76)90041-3
//' 
//'
//' @param metaW  metacommunity adyacency matrix 
//' @param m      A numeric vector of species' migration rates (probability) from the meta-web   
//' @param e      A numeric vector of species' extinction probability 
//' @param se      A numeric vector of species' secondary extinction probability 
//' @param time   Number of time steps of simulation
//' @return       A list with the final the number of species by time S, the number of links by time L, 
//'               the time series of species STime and the adjacency matrix A with effective links. 
// [[Rcpp::export]]
List metaWebNetAssemblyCT(LogicalMatrix metaW, NumericVector m, NumericVector e, NumericVector se, int time) {
  
  if( metaW.ncol()!=metaW.ncol())  
    stop("Matrix metaW have to be square");
  
  IntegerVector SL(time);      // Species with links vector, the total number of species could be bigger. 
  IntegerVector LL(time);
  
  auto rho = metaW.nrow();     // meta web should be a square matrix
  if(m.length()!=rho)
    stop("migration vector m must have the same elements that the dimensions of metaW");  
  if(e.length()!=rho)
    stop("extinction vector e must have the same elements that the dimensions of metaW");  
  
  LogicalMatrix A(rho,rho);   // Local adyacency matrix
  LogicalVector Spc(rho);     // Vector of present species
  LogicalVector Bas(rho);     // Vector of basal species
  NumericVector lambda(rho);   // Vector of scaled transition rate 
  NumericVector lambda_se(rho);   // Vector of scaled transition rate 
  NumericVector lambda_e(rho);   // Vector of scaled transition rate 
  IntegerMatrix ST(rho,time); // Vector of species populations in time
  unsigned long L=0;
  double totR = 0;
  
  // Detect Basal species
  // Build cumulative probabilities
  //
  auto r=colSums(metaW);
  for(auto s=0; s<rho; s++){
    if(r[s]==0)
      Bas[s]=true;
    totR += m[s] + e[s] + se[s];
    lambda[s] = totR;
    lambda_se[s] = lambda[s] - se[s];
    lambda_e[s] = lambda_se[s] - e[s];
  }
  lambda = lambda / totR;
  lambda_se = lambda_se / totR;
  lambda_e = lambda_e / totR;  
  
  //
  //   0 --> lambda_e[0] --> lambda_se[0] --> lambda[0] ---> lambda_e[1] --> lambda_se[1] --> lambda[1]=1
  //      |-> m[0]        |-> e[0]         |-> se[0]
  
  DBG("Basales : ",Bas)
  DBG("TotR    : ",totR)
  DBG("lambda    : ",lambda)
  DBG("lambda_se : ",lambda_se)
  DBG("lambda_e  : ",lambda_e)
    
  // Continuous time 
  //
  double T = 0;
  
  for(auto t=1; t<time; ){

    DBG("\n============ t: ",t)
    DBG("============ T: ",T)
    DBG("Spc:",Spc)
    
    
    // Add Species m*(rho-S)
    //
    auto rnd=runif(1)[0];  
    DBG("\nrandom: ",rnd)

    for(auto s=0; s<rho; ){
      if(rnd<lambda[s]){
        DBG("\nSpecie: ",s)
        
        if(rnd >= lambda_se[s]) { 
          // Secondary extinction
          DBG("\nlambda_se: ", lambda_se[s])
          
          if(!Bas[s] && Spc[s]){ 
            // if it is not basal secondary exctinction
            //
            auto colsum = 0;
            for( auto j=0; j<rho; j++) {   // If the column sum is 0 no preys 
              colsum += A(j,s);
            }
            if( colsum == 0) {  
              DBG("EXTINCION SECUNDARIA s ",s)
              for( auto j=0; j<rho; j++) {   // Delete row i (predators), column i (preys) is already 0
                  A(s,j)=false; // j-->i
              }
              Spc[s]=false;
            }
          }
        }
        else if(rnd >= lambda_e[s]) {
          DBG("\nlambda_e: ", lambda_e[s])
          // Extinction 
          // if species i extinct all column of interactions is deleted
          if( Spc[s]){
            for( auto j=0; j<rho; j++) {   
              A(j,s)=false; // j-->i
            }
            // if species i extinct all row of interactions is deleted
            for( auto j=0; j<rho; j++) {   
              A(s,j)=false; // j-->i
            }
            Spc[s]=false;             // extinction
            DBG("EXTINCION s ",s)
              
          }
        } 
        else {
          DBG("\nMigration s: ", s)
          // Migration
          //
          Spc[s]=true;             
          for( auto j=0; j<rho; j++) { 
            if( Spc[s] && Spc[j] ) {
              if( metaW(j,s) ){
                A(j,s)=true; // s<--j  Predator
              }
              if( metaW(s,j) ){
                A(s,j)=true; // s-->j  Prey
              }
            }
          }
        }
        break;
      } else s++;
    }
    
    // Update number of species
    // Update time as an exponential distribution totR

    auto deltaT = -log(runif(1)[0])/totR;
    T += deltaT;
    if(T>=t)  {
      DBG3("\nT > t ", T, " >= ", t)
      DBG("Interactions A: \n",A)
      L = sum(A);
      SL[t] = sum((colSums(A)>0)+(rowSums(A)>0)>0); // Active links
      ST(_,t) = Spc;
      
      DBG("L :", L)
      DBG("St:",SL[t])
      DBG("Spc:",Spc)
      DBG("Basal:",Bas)
      LL[t]=L;  
      t++;
    }
  }
    
  return List::create(Named("S") = SL,
                      Named("L") = LL,
                      Named("STime") = ST,
                      Named("A") = A);

}


/*** R

A <- matrix(c(c(0,1,1,1,0),c(0,0,1,0,1),c(0,0,0,0,0),c(0,0,1,0,0),c(0,0,0,1,0)),nrow = 5,byrow=TRUE)
# require(EcoNetwork)
# g <- igraph::graph_from_adjacency_matrix(A,mode="directed")
# plotTrophLevel(g,vertexLabel = TRUE, vertexSizeFactor = 10)
#
# If m=0 no species 
#
mm <- rep(1,5) #c(0.2,0.2,0.2,.2,.2)
ee <- rep(0.6, 5)
se <- rep(1, 5)
set.seed(123)

# A0 <- metaWebNetAssembly(A,mm,ee,se,10)
# A0
# A0$L[10]==1
# A0$S[10]==2

A0 <- metaWebNetAssemblyCT(A,mm,ee,se,4)
A0
A0$L[10]==1
A0$S[10]==2


# set.seed(123)
# A1 <- metaWebNetAssembly(A,mm,ee,1000)
# any(A1$L!=A0$L)

# If a=0 a full interaction matrix A is expected, upper triangle less the basal species TRUE
# A0 <- metaWebNetAssembly(A,0.2,0,20)
# A0



# If q=1 and a=0 a full interaction matrix A is expected, upper triangle less the basal species TRUE
# A <- matrix(c(c(0,1,1,1,1),c(0,0,1,1,1),c(0,0,0,1,1),c(0,0,0,0,1),c(0,1,0,0,0)),nrow = 5,byrow=TRUE)
# A0 <- metaWebNetAssembly(A,0.5,0.2,20)
# A0

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
