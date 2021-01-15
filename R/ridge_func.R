# Produces Figure 1.

library(akima)

KNN_K <- 10
scale <- 1

func <- function (x, y) {
    z <- x + y
    z <- 4*pi*(z - 0.5)
   return(sin(z) * cos(z) * exp(-z/10))
}

func2 <- function (x, y) {
    z <- sqrt(x^2+y^2)
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

pdf("images/ridge_func.pdf", width = scale*7, height = scale*7)
g <-interp(x, y, z, nx = imgrid_size, ny = imgrid_size)
image(g, col=heat.colors(128), main = 'A Ridge Function')
dev.off()

z <- rep(NA, imgrid_size**2)
for (i in 1:length(x)) {
    z[i] <- func2(x[i],y[i])
}

pdf("images/not_ridge_func.pdf", width = scale*7, height = scale*7)
g <-interp(x, y, z, nx = imgrid_size, ny = imgrid_size)
image(g, col=heat.colors(128), main = 'Not a Ridge Function')
dev.off()

