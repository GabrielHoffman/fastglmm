/***************************************************************
 * @file	fastglmm.h
 * @author	Gabriel Hoffman
 * @email	gabriel.hoffman@mssm.edu
 * @brief	Import headers
 * Copyright (C) 2024 Gabriel Hoffman
 **************************************************************/

#ifndef _FASTGLMM_H_
#define _FASTGLMM_H_

// if -D USE_R, use RcppArmadillo library
#ifdef USE_R
// [[Rcpp::depends(RcppParallel)]]  
#include <RcppArmadillo.h>
#else
#include <armadillo>
#endif

#include "fastglmm_fit.h"
#include "ModelFit.h"
#include "lmmFitResponses.h"
#include "lmmFitFeatures.h"
#include "glmmFitResponses.h"
#include "glmmFitFeatures.h"


#endif