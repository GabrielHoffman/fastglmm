

# in Rcpp pass family link functions as Function 
# as in https://github.com/jaredhuling/fastglm/blob/9a04daa4a99761fee4fc87ecdb100a530f96b161/R/fit_glm.R#L199




# Given a formula with response, random effect, and offset
# create form_fixed and form_no_offset
# formula = a ~ (1|b) + (1|c) + offset(v)
#
# process_formula(formula)
# $form_fixed
# a ~ offset(v)
#
# $form_no_offset
# a ~ (1 | b) + (1 | c)
#
#' @importFrom reformulas findbars
process_formula = function(formula, data){

	# 1) extract fixed effects
	form_fixed <- nobars(formula)
	str_rnd_only = paste0('(', as.character(findbars(formula)), ')', collapse=' + ')

	# 2) remove offset
	# adapted from MASS::glmmPQL()
	Terms <- if(missing(data)) terms(form_fixed) else terms(form_fixed, data = data)
	offt <- attr(Terms, "offset")
	offvars <- as.character(attr(Terms, "variables"))[offt + 1L]

	# if offset is the only term, make intercept explicit
  tf <- drop.terms(Terms, offt, keep.response = TRUE)
  tl <- attr(tf,"term.labels")
  if( length(tl) == 0) tl <- "1"
  formula2 = reformulate(tl, response = form_fixed[[2L]],
                      intercept = attr(tf, "intercept"),
                      env = environment(form_fixed))

  # add back random effects
  form_no_offset <- update(formula2, paste('. ~ . +', str_rnd_only))

  list( form_fixed = form_fixed, form_no_offset = form_no_offset)
}



#' Fit Generalized Linear Mixed Model via PQL
#' 
#' Fit generalized linear mixed model (GLMM) with a single random effect using penalized quasi-likelihood (PQL)
#' 
#' @param formula a two-sided linear formula object describing both the fixed-effects and random-effects part of the model, with the response on the left of a \code{~} operator and the terms, separated by \code{+} operators, on the right.  Random-effects terms are distinguished by vertical bars (\code{|}) separating expressions for design matrices from grouping factors.
#' @param data an optional data frame containing the variables named in
#' @param family a description of the error distribution and link function to be used in the model.  
#' @param weights an optional vector of prior weights with a value for each sample.  
#' @param delta if \code{NULL} estimate delta, if value is given uses this fixed value
#' @param delta.range min and max values (in log space), of the search space for delta to fit the random effect
#' @param maxit max number of PQL iterations
#' @param tol convergence criterion for the 1D search of the delta space
#' @param tol.eta convergence criterion \code{eta} in the PQL iteration
#' @param doCoxReid use Cox-Reid correction for estimating theta in negative binomial model
#' @param lambda ridge shrinkage parameter
#' @param nthreads number of threads
#'
#' @return model fit object
#' @examples
#' library(MASS)
#' library(lme4)
#' 
#' # GLMM via Laplace approximation
#' fit = glmer(y ~ trt + I(week > 2) + (1 | ID),
#'		family = binomial(), data = bacteria)
#' coef(summary(fit))
#' 
#' # GLMM via PQL
#' fit = glmmPQL(y ~ trt + I(week > 2), random = ~ 1 | ID,
#'		family = binomial, data = bacteria, verbose = FALSE)
#' coef(summary(fit))
#' 
#' # GLMM via PQL
#' fit = fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'		family = binomial(), data = bacteria)
#' coef(summary(fit))
#
#' @import stats 
#' @importFrom reformulas nobars
#' @importFrom methods is
#' @export
fastglmm = function (formula, data, family = gaussian(), weights = NULL, delta = NULL, delta.range = c(-10, 10), maxit = 100, tol = 1e-3, tol.eta = 1e-3, doCoxReid=nrow(data) < 1000, lambda = 0, nthreads = 6){

	mc <- match.call()

	# simplest way to extract data
	formula <- as.formula(formula, env = , parent.frame(1L))

	# check that formula has exactly 1 random effect
	fb <- findbars(formula)
	if (length(fb) == 0) {
		stop("formula must contain exactly 1 random effect, but none were specified")
	}
	if (length(fb) > 1) {
		stop("formula must contain exactly 1 random effect, but ", length(fb), " were specified")
	}

	# check that only 1 random effect variable is used
	vs <- all.vars(fb[[1]])
	if (length(vs) != 1) {
		stop("Only one variable can be used in the random effect")
	}

	# formula with only fixed effects
	form.fixed <- nobars(formula)

	# get variables used in response
	respVar <- all.vars(update(form.fixed, . ~ 1))
	
	## family
  if(is.character(family))
		family <- get(family, mode = "function", envir = parent.frame())
  if(is.function(family)) family <- family()
  if(is.null(family$family)) {
		print(family)
		stop("'family' not recognized")
  }
  if( maxit < 1 ){
  	stop("maxit must be >= 1")
  }

	# decompose formula
	fres <- process_formula( formula, data)
	form_fixed <- fres$form_fixed

	# drop rows with any NA values in active variables
	data_sub <- data[,colnames(data) %in% all.vars(formula)]
	data <- data[rowSums(is.na(data_sub)) == 0, colnames(data_sub)]
	data <- droplevels(data)

  # if all columns of response are in the data matrix
  # for vector response, or cbind(v1, v2)
  #   where v1 and v2 are cols in data
  # this uses the standard R processing for formulas
  if (all(respVar %in% colnames(data))) {
    # extract data
    # *very* slow when reponse is a large matrix
    mf <- model.frame(form.fixed, data, drop.unused.levels = TRUE)
    design <- model.matrix(mf, data)
    y <- model.response(mf)
    offset <- model.offset(mf)
  } else {
    # if Y is a matrix in the parent environment
    # get matrix directly from parent

    form2 <- update(form.fixed, NULL ~ .)

    # Not allowed: y^2 ~ x
    if (respVar != as.character(form.fixed)[2]) {
      stop("Function cannot be applied to reponse: ", as.character(form.fixed)[2])
    }

    mf <- model.frame(form2, data, drop.unused.levels = TRUE)
    design <- model.matrix(mf, data)
    y <- eval.parent(parse(text = respVar))
    offset <- model.offset(mf)
  }

  if( is.factor(y) ){
  	y <- as.numeric(y) - 1
  }
  if( is.character(y) ){
  	stop("Response must be numeric or factor")
  }

	if( is.null(weights) ) weights <- rep(1, nrow(data))
	if( is.null(offset) ) offset <- rep(0, nrow(data))

  Z <- preprocess_indicator(data[[vs]])
	dcmp <- indicator_decomp( Z )

	# use function based on sparse U
	if( is(dcmp$vectors, "sparseMatrix") ){
		fxn = .fastglmm_ms
	}else{
		fxn = .fastglmm_mm
	}

	# run GLMM
	fit = fxn(y = y, 
						X = design, 
						U = dcmp$vectors, 
						s = dcmp$values, 
						weights = weights, 
						offset = offset, 
						family = getFamilyString(family), 
						dcmpMethod = "categorical",
						delta = ifelse(is.null(delta), -1, delta), 
		        left = delta.range[1],
		        right = delta.range[2],
						tol = tol, 
						tol_eta = tol.eta,
						maxit = maxit,
						lambda = lambda,
						nthreads = nthreads,
						doCoxReid = doCoxReid)

	fit$s <- c(fit$s)
	fit$Z <- Z

	# format output
	fit <- as.fastlmm(fit, design = design, offset = offset, method = "PQL")

  fit$response <- y
  fit$doCoxReid <- doCoxReid

  if( grepl("^nb:", fit$family) ){
  	# convert NB string to negative.binomial(theta)
		fit$family <- stringToNbFamily( fit$family )
	}else{
		fit$family <- family
	}

	# for negative binomial
	if( isNB(family) ){
		# negative.binomial(NA) 
		# QL dispersion is 1
		if( getFamilyString(family) == "nb"){
			fit$family$dispersion <- 1
		}else{			
			# negative.binomial(x)
			# QL dispersion is estimated from data
			fit$family$dispersion <- NA
		}
	}

	fit$prior.weights <- weights
	fit$iter.pql <- fit$niter
	fit$formula <- formula	
  fit$data <- data
  fit$lambda <- lambda
	class(fit) <- c("fastglmm", "fastlmm")

	attr(fit, "call") <- mc
	fit
} 


