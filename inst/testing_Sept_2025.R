


q()
R

library(tidyverse)
library(RUnit)
library(lme4)
library(fastglmm)

implemented = c(methods(class = "fastlmm"), 
                methods(class = "fastglmm")) %>%
                gsub("\\.fast(g?)lmm", "", .) %>%
                unique %>%
                c("formula", "weights") %>%
                sort

target = methods(class = "merMod") %>%
                gsub("\\.merMod", "", .) %>%
                unique %>%
                sort

setdiff(target, implemented)

anova
# cooks.distance
# df.residual
# hatvalues
# influence
# rstudent
simulate
# deviance.residuals
# rdf
# edf
# terms

# fastlmm
##########

w = seq(nrow(sleepstudy))
w = w / mean(w)
n = nrow(sleepstudy)

fit1 <- fastlmm(Reaction ~ Days + (1 | Subject), sleepstudy, weights = w)
fit2 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy, REML = FALSE, weights = w)


H = hatvalues(fit2, TRUE)
sum(diag(H))
edf(fit1)


# does this give correct result?
# must specify test = "F"
linearHypothesis(fit1, "Days = 12")



linearHypothesis(lm(Reaction ~ Days, sleepstudy), "Days = 12")






tr = function(x) sum(diag(x))

df.residual(fit1)
n - tr(H)
n - 2*tr(H) + tr(crossprod(H))
n - 2*tr(H) + sum(H^2)


# check df.residual
checkEqualsNumeric( n - df.residual(fit1), sum(hatvalues(fit1)))

tol = 1e-5
# For each generic function
exclude = c("coef", "print", "plot", "df.residual", "extractAIC")
for(fx in setdiff(implemented, exclude)){

  cat(fx, "\n")

  # run on fastlmm fit
  res1 = get(fx)( fit1 )

  # run on lmer fit
  res2 = get(fx)( fit2 )

  if( fx %in% c("family")){
    checkEquals( res1, res2 )
  }else if( fx %in% c("formula")){
    identical(res1, res2)
  }else if( fx %in% c("terms")){
    ids = intersect(names(attributes(res1)), names(attributes(res2)))
    checkEquals(attributes(res1)[ids], attributes(res2)[ids])
  }else if( fx %in% c("ranef")){
    checkEqualsNumeric( res1[[1]], res2[[1]], tol=tol)
  }else if( fx %in% c("weights")){
    checkEqualsNumeric( c(res1), c(res2), tol=tol)
  }else if( fx %in% c("summary")){
    checkEqualsNumeric( coef(res1)[,1:3], coef(res2)[,1:3], tol=tol)
  }else{
    checkEqualsNumeric( res1, res2, tol=tol )
  }
}







w = seq(nrow(sleepstudy))
w[] = 1
w = w / mean(w)

fit1 <- fastlmm(Reaction ~ Days + (1 | Subject), sleepstudy, weights = w)
fit2 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy, REML = FALSE, weights = w)


hatvalues(fit2)
sum(hatvalues(fit2))


n = nrow(sleepstudy)
X = model.matrix( ~ Days, sleepstudy)
Z <- preprocess_indicator(sleepstudy$Subject)

Z = diag(sqrt(w)) %*% Z
X = diag(sqrt(w)) %*% X

delta = fit1$sigSq_e / fit1$sigSq_g
V = (tcrossprod(Z)/delta + diag(1, n))  
V2 = (tcrossprod(Z) + diag(delta, n)) / delta  



V_inv_X = solve(V, X)
H1 = diag(1, n) - solve(V) + V_inv_X %*% solve(crossprod(X, V_inv_X), t(V_inv_X))


V = tcrossprod(Z) + diag(delta, n)
H1 = diag(1, n) - delta * solve(V) + delta *solve(V, X) %*% solve(crossprod(X, solve(V, X)), t(solve(V, X)))






H2 = hatvalues(fit2, TRUE)

plot(H1, H2)
abline(0, 1, col="red")

sum(diag(H1))
sum(diag(H2))


dcmp = indicator_decomp(sleepstudy$Subject)
U = dcmp$vectors
U2U2t = diag(1,n) - tcrossprod(U)
U  = cbind(U, eigen(U2U2t)$vectors[,seq(n-18)])
U = as.matrix(U)
s = c(dcmp$values, rep(0, n-18))


X_til = crossprod(U, X)
X_dot = diag(1/(s+delta)) %*% X_til 
H1 = diag(1, n) - delta * U %*% diag(1/(s+delta)) %*% t(U) + delta * U %*% X_dot %*% solve(crossprod(X_til, X_dot)) %*% t(X_dot) %*% t(U)

diag(H1)

M = X_dot %*% solve(crossprod(X_til, X_dot)) %*% t(X_dot)

H1.diag = 1 - delta * (U^2) %*% (1/(s+delta)) + delta*diag(U %*% M %*% t(U))

sum(diag(H1))
sum(diag(H2))
sum(H1.diag)


# adapt to low rank Z
dcmp = indicator_decomp(sleepstudy$Subject)
U = dcmp$vectors
s = dcmp$values
V = U %*% diag(s) %*% t(U) + diag(delta, n)

diag(solve(V))

V_inv = U %*% diag(1/(s+delta)) %*% t(U) + diag(1/delta, n) - tcrossprod(U) / delta
diag(V_inv)


U^2 %*% (1/(s+delta)) + (1 - rowSums(U^2)) / delta


sum(diag(solve(V)))

sum(1/(c(dcmp$values, rep(0, n-18)) + delta))


sum(1/(dcmp$values + delta)) + (n-18) / delta





A = X / delta - U %*% diag(s/(delta*s + delta^2)) %*% crossprod(U, X)
D = solve(crossprod(A, X))

diag(A %*% D %*% t(A))



# hatvalues
# H_2 = I - V^{-1} +  V^{-1} X (X^T V^{-1} X)^{-1} X^T V^{-1}

Iu = crossprod(U,diag(1,nrow(X)))
Xu <- crossprod(U, X)
inv_s_delta <- 1 / (s + delta)
cp_X_low_I <- crossprod(X, diag(1,nrow(X))) - crossprod(Xu, Iu)
inv_s_delta_Iu <- inv_s_delta * Iu
QXI <- crossprod(Xu, inv_s_delta_Iu) + cp_X_low_I / delta
QXX <- crossprod(Xu, inv_s_delta_Xu) + cp_X_low / delta
s_expand = c(s, rep(0, nrow(X) - length(s)))
H = diag(1, nrow(X)) - diag(1/(s_expand+delta)) + X %*% solve(QXX, QXI)

sum(diag(H))

H = diag(1, nrow(X)) - solve(V) + solve(V) %*% X %*% solve(t(X) %*% solve(V, X)) %*% t(X) %*% solve(V)
sum(diag(H))





fit.lm = lm(Reaction ~ Days, sleepstudy)
V = diag(1,n)
H1 = diag(1, n) - solve(V) + solve(V, X) %*% solve(crossprod(X, solve(V, X)), t(solve(V, X)))
H2.diag = hatvalues(fit.lm)

plot(diag(H1), H2.diag)
abline(0, 1, col="red")






# delta = object$delta

# V = with(object, U %*% diag(s) %*% t(U) + diag(delta, n) )

# H = diag(1, n) - delta*solve(V) + 
#   delta * solve(V, X) %*% solve(crossprod(X,solve(V, X))) %*% t(solve(V, X))

# diag(H)
# # sum(diag(H))


# # Target: delta * diag(solve(V))
# h1 = delta*with(object, U^2 %*% (1/(s+delta)) + (1 - rowSums(U^2)) / delta)

# # Target: delta * sum(diag(solve(V)))
# with(object, delta*sum(1/(s+delta)) + (n-k)) 


# # # Target: 
# target = delta * diag(solve(V, X) %*% solve(crossprod(X,solve(V, X))) %*% t(solve(V, X)))
# # A = with(object, X / delta - (U * (s/(delta*s + delta^2))) %*% crossprod(U, X))

# A = with(object, X / delta - (U %*% diag(s/(delta*s + delta^2))) %*% crossprod(U, X))


# D = solve(crossprod(A, X))
# h2 = delta * diag(A %*% D %*% t(A))
# # h2 = delta * rowSums(A * (A %*% D))
# plot(target, h2)
# abline(0,1)


# A = with(object, X / delta - (U %*% diag(s/(delta*s + delta^2))) %*% crossprod(U, X))
# B = with(object, X / delta - U %*% ((s/(delta*s + delta^2)) * crossprod(U, X)))

# range(A-B)







