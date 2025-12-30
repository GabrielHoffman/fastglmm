
/***************************************************************
 * @file    exportToR_fastlmm.h
 * @author  Gabriel Hoffman
 * @email   gabriel.hoffman@mssm.edu
 * @brief   Export results to R as List
 * Copyright (C) 2024 Gabriel Hoffman
 **************************************************************/

#ifndef _EXPORT_TO_R_FASTLMM_
#define _EXPORT_TO_R_FASTLMM_

#ifdef USE_R

using namespace Rcpp; 
using namespace fastglmmLib;

// Depends on Rcpp::List, so define outside of class
const List toList(ModelFitLMM &res){

  List lst = List::create( 
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
                Named("iter")         = res.iter,
                Named("w.mean")       = res.w_mean,
                Named("s")            = res.s);

  // Set U as either full or sparse matrix
  if( res.isSet_U ) lst["U"] = res.U;
  if( res.isSet_Usp ) lst["U"] = res.Usp;

  // Set V as either full or sparse matrix
  if( res.isSet_V ) lst["V"] = res.V;
  if( res.isSet_Vsp ) lst["V"] = res.Vsp;

  return lst;
}

template <typename T1, typename T2, typename T3>
const List toList(fastlmm<T1, T2, T3> & fit){

  ModelFitLMM res = fit.get_result( true );

  return toList( res );
}

List toList( const vector<ModelFitLMM> &resList){
  List L = List::create();

  for(int i=0; i<resList.size(); i++){
    ModelFitLMM res = resList.at(i);
    L.push_back( toList( res ) );
  }
  return L;
}



// Depends on Rcpp::List, so define outside of class
const List toList(ModelFitGLMM &res){

  List lst = List::create( 
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
                Named("family")       = res.family,
                Named("iter")         = res.iter,
                Named("w.mean")       = res.w_mean,
                Named("s")            = res.s);

  // Set U as either full or sparse matrix
  if( res.isSet_U ) lst["U"] = res.U;
  if( res.isSet_Usp ) lst["U"] = res.Usp;

  // Set V as either full or sparse matrix
  if( res.isSet_V ) lst["V"] = res.V;
  if( res.isSet_Vsp ) lst["V"] = res.Vsp;

  return lst;
}

template <typename T1, typename T2, typename T3>
const List toList(fastglmm<T1, T2, T3> & fit){

  ModelFitGLMM res = fit.get_result();

  return toList( res );
}

#endif
#endif