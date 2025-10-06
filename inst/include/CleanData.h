/***************************************************************
 * @file    CleanData.h
 * @author  Gabriel Hoffman
 * @email   gabriel.hoffman@mssm.edu
 * @brief   Clean data with NAN values
 * Copyright (C) 2025 Gabriel Hoffman
 **************************************************************/

#ifndef CLEAN_DATA_H_
#define CLEAN_DATA_H_

#include "misc.h"

using namespace arma;
using namespace std;

namespace fastglmmLib {


/** Given Y, X and Weights for regression model, clean NAN values as follows. For rows in X with NAN values, set row values to zero.  For entries in Y with NAN values, set entries values to zero.  For Weights, set rows matching NANs in X and entries matching NAN values in Y to zero. 
 */
class CleanData {

  public:
    mat Y_clean, X_clean, W_clean;

  CleanData( const mat &Y, const mat &X, const mat &Weights )
    Y_clean(Y),
    X_clean(X),
    W_clean(Weights)  {

    // find rows in X with NAN values
    uvec idx_x = rows_with_nan(X);  
    X_clean.rows(idx_x).zeros();

    // Row rows in X with missing data, 
    // set rows in Weights to zero
    W_clean.rows(idx_x).zeros();

    // for indieces in Y with missing data
    // set indices in Weights and Y to zero
    match_NAN_zeros(Y_clean, W_clean);
    match_NAN_zeros(Y_clean, Y_clean);
  }

  mat get_Y() const {
    return Y_clean;
  }

  mat get_W() const {
    return W_clean;
  }

  mat get_X() const {
    return X_clean;
  }

  private:
  mat Y_clean, X_clean, Wsqrt;
};

}

#endif