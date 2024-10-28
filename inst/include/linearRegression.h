/***********************************************************************
 * @file		linearRegression.h
 * @author		Gabriel Hoffman
 * @email		gabriel.hoffman@mssm.edu
 * @brief		Evaluate linear regression with Armadollo library
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/


#ifndef LINEAR_REGRESSION_H_
#define LINEAR_REGRESSION_H_

struct ModelFit {
	vec coef;
	vec se;
	double df;
	string ID;
	vec hatvalues;
	ModelFit() {}
	ModelFit( const vec &coef, const vec &se, const double &df) : 
		coef(coef), se(se), df(df) {}
	ModelFit( const vec &coef, const vec &se, const double &df, const vec &hatvalues) : 
		coef(coef), se(se), df(df), hatvalues(hatvalues) {}
};

struct LMWork {
    mat Q, R, V;
    vec residuals;
    LMWork() {}
};


/** Linear regression model using QR decomposition.  Recycles QR and memory
 * 
 * @param X design matrix
 * @param y response vector
 * @param hat bool indicating if hatvalues should be returned
 * @param work LMWork workspace to store intermediate results
 * 
/ adapted from https://github.com/RcppCore/RcppArmadillo/blob/master/src/fastLm.cpp
/ https://genomicsclass.github.io/book/pages/qr_and_regression.html
*/
ModelFit lm(const arma::mat& X, const arma::colvec& y, const bool &hat = false, LMWork *work = nullptr) {

	int n = X.n_rows, k = X.n_cols;

    bool alloc_local = true;
    if( work == nullptr ){
        work = new LMWork();
    }

    // QR decomp
	qr_econ(work->Q, work->R, X);

    // back solve
	vec beta = solve(work->R, trans(work->Q) * y);

    // residuals
	work->residuals = y - X*beta;

    // std.errors of coefficients
    double rdf = n - k;
	double s2 = dot(work->residuals, work->residuals) / rdf; 

    // V = solve(t(R)*R)
    bool success = inv_sympd(work->V, trans(work->R)*work->R);

    vec stderr;
    if( success ){
    	stderr = sqrt(s2 * diagvec(work->V));
    }else{
        stderr = vec(k, fill::value(datum::nan));
        beta.fill(datum::nan);
    }
   
	ModelFit fit;
	if( hat ){
		vec hatvalues = diagvec(work->Q * trans(work->Q));
		fit = ModelFit( beta, stderr, rdf, hatvalues);
	}else{
		fit = ModelFit( beta, stderr, rdf);
	}

    if( alloc_local) delete work;

	return fit;
}


#endif