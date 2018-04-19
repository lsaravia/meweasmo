#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

//#define DBG(MSG,X) Rcpp::Rcout << MSG << X << std::endl;
#define DBG(MSG,X) 



//’ Simulation of a neutral network assembly procesess
//’
//’ @param rho metacommunity number of species
//’ @param m   migration rate (probability) from the meta-web   
//’ @param q   probability of interaction 
//’ @param a   mortality parameter for C   
//’ @param BB  number of initial species present 
//’ @return
// [[Rcpp::export]]
List neutralNetAssembly(double rho,double m, double q, double a, int time, int BB=0) {
  
  NumericVector S(time);
  IntegerVector LL(time);
  
  LogicalMatrix A(rho,rho);
  LogicalVector Spc(rho);
  LogicalVector Bas(rho);
  unsigned long L=0;

  // Basal species migration 
  if(BB==0 )
    BB=rpois(1,rho*m)[0];

  for(auto s=0; s<BB; s++){
    Spc[s]=true;
    Bas[s]=true;
    S[0]=BB;
  }      

  // Extinctions are fixed
  //
  double e = a;
  
  for(auto t=1; t<time; t++){
    DBG("\n============ t:",t)
    
    // Add Species m*(rho-S)
    //

    for(auto s=0; s<rho; s++){
      auto r=runif(1)[0];
      if(r<m)
        if(!Spc[s])
          Spc[s]=true;
    }
    S[t]=sum(Spc);
    DBG("Spc:", Spc)  
    DBG("St:",S[t])
      
    // Add interactions neutral model i=predator j=prey
    //  
    for( auto j=0; j<rho; j++) {
      for(auto i=0; i<rho; i++ ) {
        if( Spc[i] && Spc[j] && i!=j && !Bas[i] ){ // Double interactions j-->i and i-->j not allowed 
          auto r=runif(1)[0];
          if(r<q) {
            A(j,i)=true; // j-->i
          }
        }
      }
    }

    DBG("Interactions A: \n",A)
    L = sum(A);
    DBG("L ", L)
      
    // Extinctions  a L / S^2
    //
    // double e=a*L/(S[t]*S[t]);

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
          S[t]--;
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
              S[t]--;
              Spc[i]=false;
            }
          }
        }
      }
    } while (noSecondary);

    DBG("Extinctions A: \n",A)
    L = sum(A);  // Borrar porque no es necesario recalcular
    DBG("L :", L)
    DBG("St:",S[t])
    DBG("Spc:",Spc)
    DBG("Bas:",Bas)
    LL[t]=L;  
      
  }

  BB = sum(Bas*Spc);  // Basal species can become exctinct

  return List::create(Named("S") = S, 
                      Named("L") = LL,
                      Named("B") = BB,
                      Named("A") = A);
  
}



//’ Simulation of a Cascade fod web assembly procesess
//’
//’ @param rho    metacommunity number of species
//’ @param m      migration rate (probability) from the meta-web   
//’ @param q      probability of interaction 
//’ @param a      extinctions of species (not proportional to C)   
//’ @param time   time of simulation   
//’ @param BB     number of basal species if 0 will be ~poison(rho*m)   
//’ @return       A list with the final the number of species by time S, the number of links by time L, the number of basal species
//'               and the adjacency matrix A. 
// [[Rcpp::export]]
List cascadeNetAssembly(double rho,double m, double q, double a, int time, int BB=0) {
  
  NumericVector S(time);
  IntegerVector LL(time);
  
  LogicalMatrix A(rho,rho);
  LogicalVector Spc(rho);
  LogicalVector Bas(rho);
  unsigned long L=0;
  
  // Basal species migration 
  
  if(BB==0)
    BB=rpois(1,rho*m)[0];
    
  for(auto s=0; s<BB; s++){
    Spc[s]=true;
    Bas[s]=true;
    S[0]=BB;
  }      
  
  DBG("Bas",Bas)

  // Extinctions not proportional to connectance 
  // 
  double e=a;
  
  for(auto t=1; t<time; t++){
    DBG("\n============ t:",t)
    
    // Add Species m*(rho-S)
    //
    
    for(auto s=0; s<rho; s++){
      auto r=runif(1)[0];
      if(r<m)
        if(!Spc[s])
          Spc[s]=true;
    }
    // Update number of species
    S[t]=sum(Spc);
    DBG("Spc:", Spc)  
    DBG("Bas:",Bas)
    DBG("St:",S[t])
      
    // Add interactions Cascade model species i>j then i<--j i=predator j=prey
    //  
    for( auto j=0; j<rho; j++) {
      for(auto i=0; i<rho; i++ ) {
        if( Spc[i] && Spc[j] && i>j && !Bas[i]){
          auto r=runif(1)[0];
          if(r<q) {
            A(j,i)=true; // j-->i
          }
        }
      }
    }
    DBG("Interactions A: \n",A)
    L = sum(A);
    DBG("L ", L)
      
    // Extinctions  a L / S^2
    //
    //double e=a*L/(S[t]*S[t]);
    
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
          S[t]--;
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
              S[t]--;
              Spc[i]=false;
            }
          }
        }
      }
    } while (noSecondary);
    DBG("Extinctions A: \n",A)
    L = sum(A);  // Calculate L after secondary extinctions 
    DBG("L :", L)
    DBG("St:",S[t])
    DBG("Spc:",Spc)
    DBG("Bas:",Bas)
    LL[t]=L;  
    
  }
  BB = sum(Bas*Spc);  // Basal species can become exctinct
  
  return List::create(Named("S") = S, 
                      Named("L") = LL,
                      Named("B") = BB,
                      Named("A") = A);
  
}


/*** R
set.seed(123)

# Neutral Assembly model
#

# If m=0 no species
A0 <- neutralNetAssembly(20,0,1,1,20)
stopifnot(all(A0$L==0))
stopifnot(all(A0$A==0))
stopifnot(A0$B==0)


# If m=0 and a=0 forced basal species, only basal species remains without interactions
A0 <- neutralNetAssembly(20,0,1,0,20,4)
stopifnot(all(A0$L==0))
stopifnot(all(A0$A==0))
stopifnot(A0$B==4)


# If q=0 no interactions and no extinctions, only basal species survive
A0 <- neutralNetAssembly(20,1,0,0,20)
stopifnot(all(A0$L==0))
stopifnot(all(A0$A==0))
stopifnot(all(A0$B>1))


# If q=0 no interactions and no extinctions, and forced basal species, only basal species remains without interactions
A0 <- neutralNetAssembly(20,0.5,0,0,20,4)
stopifnot(all(A0$L==0))
stopifnot(all(A0$A==0))
stopifnot(all(A0$B==4))



# If q=1 and a=0 a full interaction matrix A is expected, less basal species  should be all TRUE less the diagonal
A0 <- neutralNetAssembly(20,0.5,1,0,20,4)
stopifnot(all(A0$B==4))
stopifnot(any(A0$L>0))
stopifnot(sum(A0$A) == 20*19 - 4*20 + 4)

# If m=1, q=0,  and a=1 all species migrate but go extinct including basal species
A0 <- neutralNetAssembly(20,1,0,1,20,4)
stopifnot(all(A0$L==0))
stopifnot(all(A0$A==0))
stopifnot(A0$B==0)

# If m=1, q=0,  and a=1 all species migrate but go extinct because of secondary exctinctions, except basal
A0 <- neutralNetAssembly(20,1,0,0,20,4)
stopifnot(all(A0$L==0))
stopifnot(all(A0$A==0))
stopifnot(A0$B==4)


# If m=1, q=1,and a=0 all species migrate and none go extinct
A0 <- neutralNetAssembly(20,1,1,0,20,4)
stopifnot(sum(A0$A) == 20*19 - 4*20 + 4)
stopifnot(any(A0$L==304))
stopifnot(A0$B==4)


# Cascade Assembly
#
# If m=0 no species
A0 <- cascadeNetAssembly(20,0,1,1,20)
stopifnot(all(A0$L==0))
stopifnot(all(A0$A==0))
stopifnot(A0$B==0)


# If m=0 and forced basal species, only basal species remains without interactions
A0 <- cascadeNetAssembly(20,0,1,0,20,4)
stopifnot(all(A0$L==0))
stopifnot(all(A0$A==0))
stopifnot(A0$B==4)


# If q=0 no interactions no exctinctions only basal species survive because secondary extinctions kill all the ones without preys
A0 <- cascadeNetAssembly(20,1,0,0,20)
stopifnot(all(A0$L==0))
stopifnot(all(A0$A==0))
stopifnot(all(A0$B>1))

# If q=0 no interactions and no extinctions, and forced basal species, only basal species remains without interactions
A0 <- cascadeNetAssembly(20,0.5,0,0,20,4)
stopifnot(all(A0$L==0))
stopifnot(all(A0$A==0))
stopifnot(all(A0$B==4))


# If q=1 and a=0 a full interaction matrix A is expected, upper triangle less the basal species TRUE
A0 <- cascadeNetAssembly(20,0.5,1,0,20,4)
stopifnot(all(A0$B==4))
stopifnot(any(A0$L>0))
stopifnot(sum(A0$A) == 20*19/2 - 4*3/2 )


# If m=1, q=0,  and a=1 all species migrate but go extinct including basal species
A0 <- cascadeNetAssembly(20,1,0,1,20,4)
stopifnot(all(A0$L==0))
stopifnot(all(A0$A==0))
stopifnot(all(A0$B==0))


# If m=1, q=0,  and a=1 all species migrate but go extinct because of secondary exctinctions, except basal
A0 <- cascadeNetAssembly(20,1,0,0,20,4)
stopifnot(all(A0$L==0))
stopifnot(all(A0$A==0))
stopifnot(all(A0$B==4))


# If m=1, q=1,and a=0 all species migrate and none go extinct
A0 <- cascadeNetAssembly(20,1,1,0,20,4)
stopifnot(sum(A0$A) == 20*19/2 - 4*3/2 )
stopifnot(any(A0$L==184))
stopifnot(all(A0$B==4))

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
