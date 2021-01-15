# This script produces the plots for Figure 3

library(akima)
library(hetGP)
library(activegp)

KNN_K <- 10
scale <- 1

func <- function (x, y) {
    z <- x + y
    z <- 4*pi*(z - 0.5)
   return(sin(z) * cos(z) * exp(-z/10))
}

imgrid_size <- 100
gridvals <- seq(0, 1, length= imgrid_size)
x <- rep(gridvals, rep(imgrid_size, imgrid_size))
y <- rep(gridvals, imgrid_size)
z <- rep(NA, imgrid_size**2)
for (i in 1:length(x)) {
    z[i] <- func(x[i],y[i])
}

set.seed(123)
N <- 400
P <- 2
# Doesn't seem to work for LHS?
#X <- maximinSA_LHS(matrix(runif(N*P), ncol = P))$design
X <- matrix(runif(N*P), ncol = P)
xstar <- c(0.5,0.55)
#X <- randomLHS(N,P)
Y <- sapply(1:nrow(X), function(i) func(X[i,1], X[i,2]))

dists <- sapply(1:N, function(n) sum((X[n,] - xstar)^2))
ranks <- rank(dists)
neighbors <- which(ranks<=KNN_K)

pdf("images/ridge_func.pdf", width = scale*7, height = scale*7)
g <-interp(x, y, z, nx = imgrid_size, ny = imgrid_size)
image(g, col=heat.colors(128))
dev.off()

fit <- mleHomGP(X, Y)
C <- C_GP(fit)$mat
ed <- eigen(C)
Lt <- ed$vectors %*% diag(sqrt(ed$values))

xyr <- cbind(x, y) %*% Lt
xr <- xyr[,1]
yr <- xyr[,2]
Xr <- X %*% Lt
xrstar <- xstar %*% Lt

distsr <- sapply(1:N, function(n) sum((Xr[n,] - xrstar)^2))
ranksr <- rank(distsr)
neighborsr <- which(ranksr<=KNN_K)

alpha <- 0.2 #Transparency

pdf("images/neigh_pre.pdf", width = scale*7, height = scale*7)
g <-interp(x, y, z, nx = imgrid_size, ny = imgrid_size)
image(g, col=heat.colors(128), main = "Original Input Space")
points(X[-unique(c(neighbors,neighborsr)),1], X[-unique(c(neighbors,neighborsr)),2], pch=3)
points(xstar[1], xstar[2], ,col='black', bg = 'white', lwd = 2, pch = 25, cex = 3)
points(X[neighbors,1], X[neighbors,2], col = 'black', bg = 'blue', pch = 21, cex = 3.5, lwd = 2)
transblack <- rgb(t(col2rgb('black')/255), alpha = alpha)
transgreen <- rgb(t(col2rgb('limegreen')/255), alpha = alpha)
points(X[neighborsr,1], X[neighborsr,2], col = transblack, bg = transgreen, pch = 21, cex = 2, lwd = 2)
dev.off()

#ds <- 7
#pdf("images/neigh_post.pdf", width = ds*sqrt(diff(range(xr))), height = sqrt(ds*diff(range(yr))))
pdf("images/neigh_post.pdf", width = scale*11, height = scale*7)
#pdf("images/neigh_post.pdf")
g <- interp(xr, yr, z, nx = imgrid_size, ny = imgrid_size)
image(g, col=heat.colors(128), main = "Rotated Input Space")
points(Xr[-unique(c(neighbors,neighborsr)),1], Xr[-unique(c(neighbors,neighborsr)),2], pch=3)
points(xrstar[1], xrstar[2], ,col='black', bg = 'white', lwd = 2, pch = 25, cex = 3)
points(Xr[neighborsr,1], Xr[neighborsr,2], col = 'black', bg = 'limegreen', pch = 21, cex = 3.5, lwd = 2)
transblack <- rgb(t(col2rgb('black')/255), alpha = alpha)
transblue <- rgb(t(col2rgb('blue')/255), alpha = alpha)
points(Xr[neighbors,1], Xr[neighbors,2], col = transblack, bg = transblue, pch = 21, cex = 2, lwd = 2)
dev.off()
