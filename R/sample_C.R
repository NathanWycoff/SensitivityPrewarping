# Contains functions for Kernel based estimate of sample measure C, as described in our article and Fukumizu 2014.

library(hetGP)
library(numDeriv)

#P <- 2
#N <- 1000
#sigma <- 1.2
#l <- rgamma(P,1,1)
#gpe <- 1e-8
##l <- c(0.1,1)
#
#X <- matrix(runif(N*P), ncol = P)
#y <- rowSums(X)
#a <- runif(P)
#a <- c(1,-1)
#
#K_yXya <- c(0, rep(1, N-1))
#
#k <- function(x1, x2) sigma * exp(-sum((x1-x2)^2/l))

# V[ \nalba f(a), f(X)]
# X - data locs
# a - place where gradient is evaluated.
cov_gen_dy <- function(a, X, l, sigma) {
    K_yXya <- hetGP:::cov_gen(X, matrix(a, nrow = 1), l, type = 'Gaussian')
    return(sigma * t(as.numeric(K_yXya) * t(2*(t(X) - a) / l)))
}

#xloc <- X[1,]
#grad(function(a) k(a, xloc), a)
#cov_gen_dy(a, X, l, sigma)

# V[\nabla f(a), \nabla f(a)]
cov_gen_dd <- function(a, X, l, sigma) {
    return(sigma*diag(2/l))
}

#jacobian(function(a2) grad(function(a1) k(a1, a2), a), a)
#cov_gen_dd(a, X, l, sigma)
#

# Double check nugget, especially for Vdf_df term.
C_at_a <- function(a, X, y, l, sigma, gpe, beta = 0, Ki = NULL) {
    # Generate prior variances.
    Vdf_df <- cov_gen_dd(a, X, l, sigma) 
    Vdf_f <- cov_gen_dy(a, X, l, sigma)

    if (missing(Ki)) {
        Vf_f <- sigma*hetGP::cov_gen(X, X, l, type = 'Gaussian')
        sol <- t(solve(Vf_f + diag(gpe, ncol = nrow(X), nrow = nrow(X)), t(Vdf_f)))
    } else {
        sol <- Vdf_f %*% Ki / sigma
    }

    # Get posterior moments via NCE. 
    post_mean <- sol %*% (y - beta)
    post_var <- Vdf_df - sol %*% t(Vdf_f)

    #return(post_var + tcrossprod(post_mean))
    return(post_var + tcrossprod(post_mean))
}

#xloc <- X[1,]
#eigen(C_at_a(xloc, X, y, l, sigma, gpe))$values

C_GP_empirical <- function(fit) {
    Cxs <- lapply(1:nrow(fit$X), function(n) C_at_a(fit$X[n,], fit$X, fit$Z, fit$theta, fit$nu_hat, Ki = fit$Ki))
    Cxh <- Reduce(function(x, y) x+y, Cxs) / nrow(fit$X)
}
