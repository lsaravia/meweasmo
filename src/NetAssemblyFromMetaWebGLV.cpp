#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

//#define DBG(MSG,X) Rcpp::Rcout << MSG << X << std::endl;
//#define DBG4(MSG,X,XX,Y,Z) Rcpp::Rcout << MSG << X << " " << XX << " " << Y << " " << Z << std::endl;
#define DBG(MSG,X) 
#define DBG4(MSG,X,XX,Y,Z) 



//' Simulation of an Assembly process from a Meta-web with stochastic Generalized Lotka-Volterra dynamics
//'
//' Simulates the assembly of a local community from a regional meta-web where the local population 
//' dynamics follow a stochastic generalized Lotka-Volterra (GLV) model with immigration. 
//' In contrast to [metaWebNetAssembly()] and [metaWebNetAssemblyCT()], species are not 
//' present or absent by definition: populations are integer abundances, a species arrives through 
//' immigration, and it is considered present while its abundance is greater than zero. 
//' Extinctions, including secondary extinctions of consumers that lose their resources, are not imposed 
//' by rules but emerge from the dynamics (demographic stochasticity and interactions).
//' 
//' The abundance \eqn{N_i} of species \eqn{i} changes in time according to
//' \deqn{\frac{dN_i}{dt} = N_i \left( r_i + \sum_j A_{ij} N_j \right) + m_i}{dN_i/dt = N_i (r_i + sum_j A_ij N_j) + m_i}
//' where \eqn{r_i} is the intrinsic growth rate (`r`), \eqn{A_{ij}}{A_ij} the interaction matrix (`metaW`) 
//' and \eqn{m_i} a constant immigration rate (`m`). 
//' 
//' The model is integrated with the tau-leap method \[1\]: time is divided in steps of size `tau`, 
//' and in each step the abundance of each species changes by a Poisson distributed number of individuals 
//' with mean `abs(rate) * tau`, where `rate` is the right hand side of the equation above; 
//' the change is positive if `rate > 0` and negative otherwise. Abundances are truncated at zero. 
//' Species are updated sequentially within a step, so a species uses the already updated abundances 
//' of the species before it. 
//' 
//' The signs of the interaction matrix define the type of interaction between species `i` and `j`:
//' * `metaW[i, j] < 0` and `metaW[j, i] > 0`: predator-prey, `j` is the predator and `i` the prey.
//' * `metaW[i, j] > 0` and `metaW[j, i] > 0`: mutualism.
//' * `metaW[i, j] < 0` and `metaW[j, i] < 0`: competition.
//' * `metaW[i, j] != 0` and `metaW[j, i] == 0`: commensalism (if positive) or amensalism (if negative).
//' * `metaW[i, i] < 0`: self-limitation (density dependence). Without it, or without negative 
//'   interactions, abundances can grow without bound.
//' 
//' Parameters of this kind of matrix can be generated at random with [generateRandomGLVadjMat()] 
//' (structure) and [generateGLVparmsFromAdj()] (interaction strengths, `r` and `m`). 
//' 
//' @references
//' 
//' 1. Gillespie, D. T. (2001). Approximate accelerated stochastic simulation of chemically reacting systems. 
//' Journal of Chemical Physics, 115(4), 1716–1733. https://doi.org/10.1063/1.1378322
//'
//' @param metaW  Numeric square matrix with the interaction coefficients \eqn{A_{ij}}{A_ij} of the meta-web: 
//'               the effect of species `j` on the per capita growth rate of species `i`, in 
//'               `metaW[i, j]`. The values and the signs are the parameters of the GLV model, 
//'               the diagonal is the self-limitation.
//' @param m      Numeric vector of length `nrow(metaW)`, immigration rate of each species, in individuals 
//'               per unit of time. It is not a probability: a species with `m > 0` receives immigrants 
//'               continuously, even when it is locally absent. Use 0 for species that should not migrate.
//' @param r      Numeric vector of length `nrow(metaW)`, intrinsic per capita growth rate of each species. 
//'               It could be negative (e.g. consumers that need resources to persist).
//' @param ini    Numeric vector with the initial abundance of each species (length `nrow(metaW)`), or a single
//'               value that is recycled to all species (e.g. 0 for an initially empty community).
//'               Values are truncated to integers.
//' @param time   Number of time points of the simulation. The initial condition is time point 1 
//'               and the state is recorded at each time unit, so the simulation covers `time - 1` time units.
//' @param tau    Size of the integration step of the tau-leap method, in time units. Each time unit is divided 
//'               in `1/tau` steps, so `1/tau` should be an integer. Smaller values are more accurate but slower.
//' @return       A list with:
//' * `S`: integer vector of length `time`, number of species present (abundance > 0) at each time point.
//' * `L`: integer vector of length `time`, number of directed links at each time point, that is 
//'   the number of non-zero interactions of `metaW` between present species; a predator-prey 
//'   interaction counts as one link (from prey to predator), other interactions count in each direction 
//'   in which they are non-zero.
//' * `STime`: integer matrix of dimension `nrow(metaW)` x `time` with the abundance of each species 
//'   (rows) at each time point (columns).
//' * `A`: logical adjacency matrix of the links in the last time point, `A[i, j]` is `TRUE` if 
//'   there is a link from `i` to `j` (for predator-prey from the prey to the predator).
//' @seealso [metaWebNetAssembly()] and [metaWebNetAssemblyCT()] for models with presence/absence 
//' dynamics, [generateGLVparmsFromAdj()] to generate the parameters.
//' @examples
//' # Predator-prey dynamics with immigration
//' A <- matrix(c(-0.001, -0.08,
//'                0.01,   0   ), nrow = 2, byrow = TRUE)
//' set.seed(123)
//' out <- metaWebNetAssemblyGLV(A, m = c(0.1, 0.2), r = c(1.5, -0.25), ini = c(0, 0), 
//'                              time = 50, tau = 0.1)
//' matplot(t(out$STime), type = "l", lty = 1, xlab = "Time", ylab = "Abundance")
//' out$S
//' 
//' # Random meta-web of 20 species, starting from an empty community
//' set.seed(1)
//' adj <- generateRandomGLVadjMat(20, 0.2, c(0.1, 0.3, 0.3, 0.15, 0.15))
//' parms <- generateGLVparmsFromAdj(adj, ef = 0.1, predIntAvg = 0.01, migrAvg = 0.01)
//' out <- metaWebNetAssemblyGLV(parms$interM, parms$m, parms$r, ini = 0, time = 100, tau = 0.1)
//' plot(out$S, type = "l", xlab = "Time", ylab = "Number of species")
//' @export
// [[Rcpp::export]]
List metaWebNetAssemblyGLV(NumericMatrix metaW,NumericVector m, NumericVector r, NumericVector ini, int time, double tau=0.01) {
  
  IntegerVector SL(time);      // Species with links vector, the total number of species could be bigger. 

  IntegerVector LL(time);
  auto rho = metaW.nrow();     // meta web should be a square matrix
  
  LogicalMatrix A(rho,rho);   // Local adyacency matrix
  IntegerVector Spc(rho);     // Vector of species populations
  IntegerMatrix ST(rho,time); // Vector of species populations in time

  // Initial values: a vector of length nrow(metaW), or a single value recycled to all species
  if(ini.size()==1) {
    std::fill(Spc.begin(), Spc.end(), (int) ini[0]);
  } else if(ini.size()==rho) {
    for(auto i=0; i<rho; i++) Spc[i] = (int) ini[i];
  } else {
    stop("ini must have length 1 or nrow(metaW)");
  }
  if(m.size()!=rho || r.size()!=rho)
    stop("m and r must have length nrow(metaW)");

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
