/***********************************************************************
 * @file    exportToR.h
 * @author    Gabriel Hoffman
 * @email   gabriel.hoffman@mssm.edu
 * @brief   Functions to export results to R
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/


#ifndef EXPORT_TO_R_H_
#define EXPORT_TO_R_H_

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace Rcpp;
using namespace arma;
using namespace fastlmmLib;

// List toList( const vector<ModelFit> & fitList){

//     int ncoef = fitList[0].coef.n_elem;
//     int nmodels = fitList.size(); 

//     mat coefMat(nmodels,ncoef);
//     mat seMat(nmodels,ncoef);
//     NumericVector rdf(nmodels);
//     NumericVector sigSq(nmodels);
//     vector<string> ID;
//     ID.reserve(nmodels);

//     // check optional properties of model
//     bool useVCOV = (fitList[0].vcov.n_rows !=0);
//     bool useResid = (fitList[0].residuals.n_elem !=0);
//     bool useHat = (fitList[0].hatvalues.n_elem !=0);

//     mat vcovStacked;
//     mat residuals, hatvalues;

//     if( useVCOV ){
//       // column i stores the vcov matrix for response i
//       vcovStacked = mat(pow(ncoef,2), nmodels);
//     }

//     if( useResid ){
//       residuals = mat(fitList[0].residuals.n_elem, nmodels);
//     }

//     if( useHat ){
//       hatvalues = mat(fitList[0].hatvalues.n_elem, nmodels);
//     }

//     // for each model
//     for(int i=0; i< fitList.size(); i++){
//       coefMat.row(i) = fitList[i].coef.t();
//       seMat.row(i) = fitList[i].se.t();
//       sigSq(i) = fitList[i].sigSq;
//       rdf(i) = fitList[i].rdf;
//       ID.push_back(fitList[i].ID);
//       if( useVCOV ) vcovStacked.col(i) = fitList[i].vcov.as_col();
//       if( useResid ) residuals.col(i) = fitList[i].residuals;
//       if( useHat ) hatvalues.col(i) = fitList[i].hatvalues;
//     }

//     // Create standard return values
//     CharacterVector IDcv(ID.begin(), ID.end());

//     sigSq.attr("names") = IDcv;
//     rdf.attr("names") = IDcv;

//     List lst = List::create(
//         Named("coef") = coefMat,
//         Named("se") = seMat,
//         Named("sigSq") = sigSq,
//         Named("rdf") = rdf
//       );

//     // Insert other return values
//     if( useVCOV ) lst["vcovStacked"] = vcovStacked;
//     if( useResid ) lst["residuals"] = residuals;
//     if( useHat ) lst["hatvalues"] = hatvalues;

//     // Set names
//     rownames(lst["coef"]) = IDcv;
//     rownames(lst["se"]) = IDcv;;
//     if( useVCOV )   colnames(lst["vcovStacked"]) = IDcv;
//     if( useResid )  colnames(lst["residuals"]) = IDcv;
//     if( useHat )    colnames(lst["hatvalues"]) = IDcv;

//     return lst;
// }

#endif