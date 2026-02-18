

# https://chatgpt.com/share/698f66b1-88b0-800b-9001-e28f444f6015
est_hessian = function(fit){

  delta = fit$delta
  sigSq_g = fit$sigSq_g
  s = fit$s

  n = length(fit$y)
  r = length(s)
  H = matrix(0, 2,2)

  A = sum(1/(s + delta)) + (n-r)/delta
  B = sum(1/(s + delta)^2) + (n-r)/delta^2

  H[1,1] = (n - 2*delta*A + delta^2*B) / (2*sigSq_g^2)
  H[1,2] = H[2,1] = (A - delta*B) / (2*sigSq_g^2)
  H[2,2] = B / (2*sigSq_g^2)
  H
}

est_gradient = function(fit, L){
  X = fit$design
  W = with(fit, solve(tcrossprod(Z) + diag(delta, nrow(Z))))

  # evalaute A and B faster with decorrelate identities
  A = crossprod(X, W) %*% X
  B = crossprod(X, W %*% W) %*% X

  lapply(seq(nrow(L)), function(i){

    g = c(0, 0)
    invAL = solve(A, L[i,])
    C = invAL %*% B %*% invAL
    g[1] = crossprod(L[i,], invAL) - fit$delta*C
    g[2] = C
    g
  })
}


ddf = function(fit, L = diag(1, length(coef(fit)))){

  hess = est_hessian(fit)
  A = solve(hess)
  g = est_gradient(fit, L)

  sapply(seq(nrow(L)), function(i){
    var_Lbeta = crossprod(L[i,], vcov(fit)) %*% L[i,]
    v_numerator <- 2 * var_Lbeta^2
    v_denom = crossprod(g[[i]], A) %*% g[[i]]

    v_numerator / v_denom  
  })
}

