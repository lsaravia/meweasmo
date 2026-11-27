#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

//' Simulation of an assembly process from a multilayer meta-web: C++ core
//'
//' Low level simulation used by [metaWebNetAssemblyML()]. Most users should call [metaWebNetAssemblyML()], which
//' takes the network as layers and the effects as the parameters `beta`; this function takes instead the three
//' matrices built from them by [build_multilayer_assembly()]. Use it directly to run many simulations of the same
//' network without rebuilding the matrices.
//'
//' It is the continuous time model of [metaWebNetAssemblyCT()] with two additions:
//'
//' * **Non-trophic modifiers of extinction**: the extinction rate of species *i* depends on which species
//'   are present, \eqn{e_i(t) = e_i \exp(\sum_j M_{ij} x_j(t))}, where \eqn{x_j(t)} is 1 if species *j* is
//'   present. `M` is in units of log extinction rate: [build_multilayer_assembly()] puts \eqn{-\beta_+} (divided by
//'   the number of positive partners if `normalize = TRUE`) where *j* has a positive effect on *i*, and
//'   \eqn{+\beta_-} where it has a negative one. So a negative \eqn{M_{ij}} (a mutualist, a facilitator) lowers the
//'   extinction rate of *i* while *j* is present and a positive one (a competitor) raises it.
//'   With `M = 0` the model is [metaWebNetAssemblyCT()].
//' * **Obligate dependence**: if species *i* depends on at least one species (`Dep[i, ]` has a `TRUE`) and none of
//'   them is present, *i* goes extinct with the secondary extinction rate *se*, as a predator without prey does.
//'
//' The state dependent extinction rate is simulated exactly by thinning: extinction events are proposed with the
//' upper bound \eqn{e_i \exp(\sum_j \max(M_{ij},0))} and accepted with probability \eqn{e_i(t)} divided by that
//' bound. Large positive values in `M` make the bound large and the simulation slower (more rejected proposals).
//'
//' @param metaW  Logical trophic adjacency matrix of the meta-web, `metaW[j, i]` is `TRUE` if *j* is a prey of *i*.
//' @param M      Numeric matrix, `M[i, j]` is the effect of the presence of *j* on the log extinction rate of *i*
//'   (the diagonal is ignored).
//' @param Dep    Logical matrix, `Dep[i, j]` is `TRUE` if *i* depends obligately on *j*.
//' @param m      A numeric vector of species' migration rates from the meta-web
//' @param e      A numeric vector of species' baseline extinction rates
//' @param se     A numeric vector of species' secondary extinction rates (no prey, or no obligate partner)
//' @param time   Simulated time, the state is recorded at integer times 1 to time-1
//' @return A list as [metaWebNetAssemblyCT()]: `S` and `L` (number of species with trophic links and of trophic
//'   links at each time), `STime` (species presence, species x time) and `A` (trophic adjacency matrix of the final
//'   local web), plus `rejected`, the number of rejected extinction proposals.
//' @seealso [metaWebNetAssemblyML()], [build_multilayer_assembly()], [metaWebNetAssemblyCT()]
// [[Rcpp::export]]
List metaWebNetAssemblyMLCT(LogicalMatrix metaW, NumericMatrix M, LogicalMatrix Dep,
                            NumericVector m, NumericVector e, NumericVector se, int time) {

  auto rho = metaW.nrow();
  if(metaW.ncol()!=rho) stop("Matrix metaW has to be square");
  if(M.nrow()!=rho || M.ncol()!=rho) stop("M must have the same dimensions as metaW");
  if(Dep.nrow()!=rho || Dep.ncol()!=rho) stop("Dep must have the same dimensions as metaW");
  if(m.length()!=rho || e.length()!=rho || se.length()!=rho)
    stop("m, e and se must have one element per species of metaW");

  IntegerVector SL(time);
  IntegerVector LL(time);
  LogicalMatrix A(rho,rho);
  LogicalVector Spc(rho);
  LogicalVector Bas(rho);
  LogicalVector HasDep(rho);
  NumericVector emax(rho);
  NumericVector lambda(rho), lambda_se(rho), lambda_e(rho);
  IntegerMatrix ST(rho,time);
  double rejected = 0;

  auto r = colSums(metaW);
  double totR = 0;
  for(auto s=0; s<rho; s++){
    Bas[s] = (r[s]==0);
    double up = 0;
    for(auto j=0; j<rho; j++){
      if(j!=s && M(s,j)>0) up += M(s,j);
      if(Dep(s,j)) HasDep[s] = true;
    }
    emax[s] = e[s]*std::exp(up);
    totR += m[s] + emax[s] + se[s];
    lambda[s] = totR;
    lambda_se[s] = lambda[s] - se[s];
    lambda_e[s] = lambda_se[s] - emax[s];
  }
  lambda = lambda / totR;
  lambda_se = lambda_se / totR;
  lambda_e = lambda_e / totR;

  double T = 0;
  for(auto t=1; t<time; ){
    auto rnd = runif(1)[0];
    for(auto s=0; s<rho; ){
      if(rnd<lambda[s]){
        if(rnd >= lambda_se[s]) {
          // Secondary extinction: no prey (non-basal) or no obligate partner
          if(Spc[s]){
            bool noPrey = false;
            if(!Bas[s]){
              int colsum = 0;
              for(auto j=0; j<rho; j++) colsum += A(j,s);
              noPrey = (colsum==0);
            }
            bool noPartner = false;
            if(HasDep[s]){
              noPartner = true;
              for(auto j=0; j<rho; j++) if(Dep(s,j) && Spc[j]) { noPartner = false; break; }
            }
            if(noPrey || noPartner){
              for(auto j=0; j<rho; j++) { A(s,j)=false; A(j,s)=false; }
              Spc[s]=false;
            }
          }
        }
        else if(rnd >= lambda_e[s]) {
          // Extinction with state dependent rate, by thinning
          if(Spc[s]){
            double lg = 0;
            for(auto j=0; j<rho; j++) if(j!=s && Spc[j]) lg += M(s,j);
            double pacc = (emax[s]>0) ? e[s]*std::exp(lg)/emax[s] : 0;
            if(runif(1)[0] < pacc){
              for(auto j=0; j<rho; j++) { A(j,s)=false; A(s,j)=false; }
              Spc[s]=false;
            } else rejected++;
          }
        }
        else {
          // Migration
          Spc[s]=true;
          for(auto j=0; j<rho; j++) {
            if(Spc[j]) {
              if(metaW(j,s)) A(j,s)=true;
              if(metaW(s,j)) A(s,j)=true;
            }
          }
        }
        break;
      } else s++;
    }

    T += -std::log(runif(1)[0])/totR;
    while(T>=t && t<time) {
      LL[t] = sum(A);
      SL[t] = sum((colSums(A)>0)+(rowSums(A)>0)>0);
      ST(_,t) = Spc;
      t++;
    }
  }

  return List::create(Named("S") = SL,
                      Named("L") = LL,
                      Named("STime") = ST,
                      Named("A") = A,
                      Named("rejected") = rejected);
}
