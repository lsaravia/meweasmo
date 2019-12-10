// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// neutralNetAssembly
List neutralNetAssembly(double rho, double m, double q, double a, int time, int BB);
RcppExport SEXP _meweasmo_neutralNetAssembly(SEXP rhoSEXP, SEXP mSEXP, SEXP qSEXP, SEXP aSEXP, SEXP timeSEXP, SEXP BBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type time(timeSEXP);
    Rcpp::traits::input_parameter< int >::type BB(BBSEXP);
    rcpp_result_gen = Rcpp::wrap(neutralNetAssembly(rho, m, q, a, time, BB));
    return rcpp_result_gen;
END_RCPP
}
// cascadeNetAssembly
List cascadeNetAssembly(double rho, double m, double q, double a, int time, int BB);
RcppExport SEXP _meweasmo_cascadeNetAssembly(SEXP rhoSEXP, SEXP mSEXP, SEXP qSEXP, SEXP aSEXP, SEXP timeSEXP, SEXP BBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type time(timeSEXP);
    Rcpp::traits::input_parameter< int >::type BB(BBSEXP);
    rcpp_result_gen = Rcpp::wrap(cascadeNetAssembly(rho, m, q, a, time, BB));
    return rcpp_result_gen;
END_RCPP
}
// metaWebNetAssembly
List metaWebNetAssembly(LogicalMatrix metaW, NumericVector m, NumericVector e, int time);
RcppExport SEXP _meweasmo_metaWebNetAssembly(SEXP metaWSEXP, SEXP mSEXP, SEXP eSEXP, SEXP timeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< LogicalMatrix >::type metaW(metaWSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type m(mSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type e(eSEXP);
    Rcpp::traits::input_parameter< int >::type time(timeSEXP);
    rcpp_result_gen = Rcpp::wrap(metaWebNetAssembly(metaW, m, e, time));
    return rcpp_result_gen;
END_RCPP
}
// metaWebNetAssemblyGLV
List metaWebNetAssemblyGLV(NumericMatrix metaW, NumericVector m, NumericVector r, NumericVector ini, int time, double tau);
RcppExport SEXP _meweasmo_metaWebNetAssemblyGLV(SEXP metaWSEXP, SEXP mSEXP, SEXP rSEXP, SEXP iniSEXP, SEXP timeSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type metaW(metaWSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type m(mSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ini(iniSEXP);
    Rcpp::traits::input_parameter< int >::type time(timeSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(metaWebNetAssemblyGLV(metaW, m, r, ini, time, tau));
    return rcpp_result_gen;
END_RCPP
}
// calcPropInteractionsGLVadjMat
NumericVector calcPropInteractionsGLVadjMat(NumericMatrix adjM, NumericVector spc);
RcppExport SEXP _meweasmo_calcPropInteractionsGLVadjMat(SEXP adjMSEXP, SEXP spcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type adjM(adjMSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type spc(spcSEXP);
    rcpp_result_gen = Rcpp::wrap(calcPropInteractionsGLVadjMat(adjM, spc));
    return rcpp_result_gen;
END_RCPP
}
// generateGLVparmsFromAdj
List generateGLVparmsFromAdj(NumericMatrix adjM, double ef, double predIntAvg, NumericVector selfLimAvg, double migrAvg, int preserveInt);
RcppExport SEXP _meweasmo_generateGLVparmsFromAdj(SEXP adjMSEXP, SEXP efSEXP, SEXP predIntAvgSEXP, SEXP selfLimAvgSEXP, SEXP migrAvgSEXP, SEXP preserveIntSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type adjM(adjMSEXP);
    Rcpp::traits::input_parameter< double >::type ef(efSEXP);
    Rcpp::traits::input_parameter< double >::type predIntAvg(predIntAvgSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type selfLimAvg(selfLimAvgSEXP);
    Rcpp::traits::input_parameter< double >::type migrAvg(migrAvgSEXP);
    Rcpp::traits::input_parameter< int >::type preserveInt(preserveIntSEXP);
    rcpp_result_gen = Rcpp::wrap(generateGLVparmsFromAdj(adjM, ef, predIntAvg, selfLimAvg, migrAvg, preserveInt));
    return rcpp_result_gen;
END_RCPP
}
// generateRandomGLVadjMat
IntegerMatrix generateRandomGLVadjMat(int numSp, double C, NumericVector propInt);
RcppExport SEXP _meweasmo_generateRandomGLVadjMat(SEXP numSpSEXP, SEXP CSEXP, SEXP propIntSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type numSp(numSpSEXP);
    Rcpp::traits::input_parameter< double >::type C(CSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type propInt(propIntSEXP);
    rcpp_result_gen = Rcpp::wrap(generateRandomGLVadjMat(numSp, C, propInt));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_meweasmo_neutralNetAssembly", (DL_FUNC) &_meweasmo_neutralNetAssembly, 6},
    {"_meweasmo_cascadeNetAssembly", (DL_FUNC) &_meweasmo_cascadeNetAssembly, 6},
    {"_meweasmo_metaWebNetAssembly", (DL_FUNC) &_meweasmo_metaWebNetAssembly, 4},
    {"_meweasmo_metaWebNetAssemblyGLV", (DL_FUNC) &_meweasmo_metaWebNetAssemblyGLV, 6},
    {"_meweasmo_calcPropInteractionsGLVadjMat", (DL_FUNC) &_meweasmo_calcPropInteractionsGLVadjMat, 2},
    {"_meweasmo_generateGLVparmsFromAdj", (DL_FUNC) &_meweasmo_generateGLVparmsFromAdj, 6},
    {"_meweasmo_generateRandomGLVadjMat", (DL_FUNC) &_meweasmo_generateRandomGLVadjMat, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_meweasmo(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
