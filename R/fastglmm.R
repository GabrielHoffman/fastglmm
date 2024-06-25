

fastglmmPQL = function (fixed, random, family, data, correlation, weights, control, niter = 10, verbose = TRUE, ...) {
    if (!requireNamespace("nlme", quietly = TRUE)) 
        stop("package 'nlme' is essential")
    if (is.character(family)) 
        family <- get(family)
    if (is.function(family)) 
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    m <- mcall <- Call <- match.call()
    nm <- names(m)[-1L]
    keep <- is.element(nm, c("weights", "data", "subset", "na.action"))
    for (i in nm[!keep]) m[[i]] <- NULL
    allvars <- if (is.list(random)) 
        c(all.vars(fixed), names(random), unlist(lapply(random, 
            function(x) all.vars(formula(x)))))
    else c(all.vars(fixed), all.vars(random))
    Terms <- if (missing(data)) 
        terms(fixed)
    else terms(fixed, data = data)
    offt <- attr(Terms, "offset")
    have_offset <- length(offt) >= 1L
    if (have_offset) {
        offvars <- as.character(attr(Terms, "variables"))[offt + 
            1L]
        allvars <- c(allvars, offvars)
    }
    if (!missing(correlation) && !is.null(attr(correlation, "formula"))) 
        allvars <- c(allvars, all.vars(attr(correlation, "formula")))
    Call$fixed <- eval(fixed)
    Call$random <- eval(random)
    m$formula <- as.formula(paste("~", paste(allvars, collapse = "+")))
    environment(m$formula) <- environment(fixed)
    m$drop.unused.levels <- TRUE
    m[[1L]] <- quote(stats::model.frame)
    mf <- eval.parent(m)
    off <- model.offset(mf)
    if (is.null(off)) 
        off <- 0
    fixed2 <- if (have_offset) {
        tf <- drop.terms(Terms, offt, keep.response = TRUE)
        reformulate(attr(tf, "term.labels"), response = fixed[[2L]], 
            intercept = attr(tf, "intercept"), env = environment(fixed))
    }
    else fixed
    wts <- model.weights(mf)
    if (is.null(wts)) 
        wts <- rep(1, nrow(mf))
    mf$wts <- wts
    fit0 <- glm(formula = fixed, family = family, data = mf, 
        weights = wts, ...)
    w <- fit0$prior.weights
    eta <- fit0$linear.predictors
    zz <- eta + fit0$residuals - off
    wz <- fit0$weights
    fam <- family
    nm <- names(mcall)[-1L]
    keep <- is.element(nm, c("fixed", "random", "data", "subset","na.action", "control"))
    for (i in nm[!keep]) mcall[[i]] <- NULL
    fixed2[[2L]] <- quote(zz)
    mcall[["fixed"]] <- fixed2
    mcall[[1L]] <- quote(nlme::lme)
    mcall$random <- random
    mcall$method <- "ML"
    if (!missing(correlation)) 
        mcall$correlation <- correlation
    mcall$weights <- quote(nlme::varFixed(~invwt))
    mf$zz <- zz
    mf$invwt <- 1/wz
    mcall$data <- mf
    for (i in seq_len(niter)) {
        if (verbose) 
            message(gettextf("iteration %d", i), domain = NA)
        fit <- eval(mcall)
        etaold <- eta
        eta <- fitted(fit) + off
        if (sum((eta - etaold)^2) < 1e-06 * sum(eta^2)) 
            break
        mu <- fam$linkinv(eta)
        mu.eta.val <- fam$mu.eta(eta)
        mf$zz <- eta + (fit0$y - mu)/mu.eta.val - off
        wz <- w * mu.eta.val^2/fam$variance(mu)
        mf$invwt <- 1/wz
        mcall$data <- mf
    }
    attributes(fit$logLik) <- NULL
    fit$call <- Call
    fit$family <- family
    fit$logLik <- as.numeric(NA)
    if (have_offset) {
        fit$fitted <- fit$fitted + off
        fit$have_offset <- TRUE
        fit$fixed2 <- fixed2
        fit$offvars <- offvars
    }
    oldClass(fit) <- c("glmmPQL", oldClass(fit))
    fit
}
