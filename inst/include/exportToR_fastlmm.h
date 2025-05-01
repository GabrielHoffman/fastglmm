
/***************************************************************
 * @file    exportToR_fastlmm.h
 * @author  Gabriel Hoffman
 * @email   gabriel.hoffman@mssm.edu
 * @brief   Export results to R as List
 * Copyright (C) 2024 Gabriel Hoffman
 **************************************************************/

#ifndef _EXPORT_TO_R_FASTLMM_
#define _EXPORT_TO_R_FASTLMM_

// #include "fastglmm.h"

using namespace Rcpp; 
using namespace fastglmmLib;

// Depends on Rcpp::List, so define outside of class
const List toList(ModelFitLMM &res){

  return List::create( 
                Named("logLik")       = res.logLik, 
                Named("coefficients") = res.coef,
                Named("se")           = res.se,
                Named("vcov")         = res.vcov, 
                Named("weights")      = res.weights,
                Named("delta")        = res.delta,
                Named("sigSq_g")      = res.sigSq_g,
                Named("sigSq_e")      = res.sigSq_e,
                Named("ru")           = res.ru,
                Named("y")            = res.y,
                Named("iter")         = res.iter);
}

template <typename T1, typename T2, typename T3>
const List toList(fastlmm<T1, T2, T3> & fit){

  ModelFitLMM a = fit.get_result();

  return toList( a);
}


List toList( const vector<ModelFitLMM> &resList){
  List L = List::create();

  for(int i=0; i<resList.size(); i++){
    ModelFitLMM a = resList.at(i);
    L.push_back( toList(a) );
  }
  return L;
}

#endif