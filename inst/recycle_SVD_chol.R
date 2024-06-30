


library(Rfast)
library(Matrix)

n = 2000
p = 4000
X = matrnorm(n,p)


dcmp = svd(X, nv=0) # UDV^T


D = Diagonal(length(dcmp$d),dcmp$d)
W = Diagonal(n, seq(n))

system.time(
	dcmp2 <- svd(W %*% X)
)

W %*% X
W %*% dcmp$u %*% diag(dcmp$d) %*% t(dcmp$v)

system.time(
res <- qr(W %*% dcmp$u)
)

# can do this implicitly
A = qr.Q(res) %*% qr.R(res) %*% diag(dcmp$d) %*% t(dcmp$v)

W %*% X
A

system.time(
	dcmp3 <- svd( W %*% dcmp$u %*% D, nv=0 )
)

dcmp3$u
svd(W %*% X)$u

dcmp3$d
svd(W %*% X)$d



system.time({
	K <- crossprod(W %*% dcmp$u %*% D)
	edecomp <- eigen(K, symmetric=TRUE)
})

r = runif(ncol(dcmp$u))
system.time({
ch = Cholesky(crossprod(W %*% dcmp$u %*% D))
})

system.time({
	replicate(100, 
	res <- solve(ch, r, Imult=32))
})

A = W %*% dcmp$u %*% D
A = as(A, "sparseMatrix")
XX = crossprod(A)
ch1 = Cholesky(XX)
u1 <- update(ch1, XX, mult=1)
u2 <- update(ch1, t(A), mult=1) 
stopifnot(all.equal(u1,u2, tol=1e-14))

A = dcmp$u %*% D
A = as(A, "sparseMatrix")
XX = crossprod(A)
ch1 = Cholesky(XX)
u3 <- update(ch1, crossprod(W %*% A), mult=1)
u4 <- update(ch1, t(W %*% A), mult=1) 
stopifnot(all.equal(u3,u4, tol=1e-14))
stopifnot(all.equal(u1,u4, tol=1e-14))

system.time(u4 <- update(ch1, t(W %*% A), mult=1) )
G = t(W %*% A)
system.time(u4 <- update(ch1, G, mult=1) )


# https://www.rdocumentation.org/packages/Matrix/versions/1.4-0/topics/CHMfactor-class

u1 <- update(CX, XX,    mult=pi)


u2 <- update(CX, t(M1), mult=pi) 
stopifnot(all.equal(u1,u2, tol=1e-14))

# Fast update of Cholesky with weighting
D = Diagonal(nrow(M1), seq(nrow(M1)))
u1 <- update(CX, crossprod(D %*% M1),    mult=pi)
u2 <- update(CX, t(D %*% M1), mult=pi) 
stopifnot(all.equal(u1,u2, tol=1e-14))








