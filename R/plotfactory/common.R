#!/usr/bin/Rscript
#  R/plotfactory/common.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 12.16.2020

cex_axis <- 1.4

#agp_inds <- c(1, 4, 10, 16, 5, 11, 17)
agp_inds <- grep("aGP", method_names)
vecc_inds <- grep("Vecchia", method_names, ignore.case = T)
knn_inds <- grep("knn", method_names)
cols <- c("white", rep('maroon', length(agp_inds)), rep('lightblue', length(vecc_inds)), rep('violet', length(knn_inds)))
slightly_darker_cols <- rgb(t(col2rgb(cols)/255)*0.6)

oma <- c(0, 0, 1.4, 0)
mar <- c(5, 4, 0.5, 0.4) + 0.1
cex_axis <- 1.4

names <- gsub('lengthscale','B', method_names)
names <- gsub('lebesgue','L', names)
names <- gsub('empirical','S', names)
names <- gsub('trunc','T', names)
names <- gsub('_','-', names)
names <- gsub('subset','s', names)
names <- gsub('aGP','laGP', names)
names <- gsub('knn','KNN', names)
names <- gsub('Vecchia','vecc', names, ignore.case = T)
#cols <- rep(c("white", "maroon", "lightblue", "violet"), rep(2,4))

height_scale <- 1.3

reorder <- c(2, agp_inds, vecc_inds, knn_inds)

mses <- mses[,,reorder,drop=F]
scores <- scores[,,reorder,drop=F]
names <- names[reorder]

notna <- !apply(mses, 3, function(x) any(is.na(x)))

#mses <- mses[,,notna,drop=F]
#scores <- scores[,,notna,drop=F]
#names <- names[notna]
#cols <- cols[notna]
#slightly_darker_cols <- slightly_darker_cols[notna]

hasscore <- setdiff(1:length(names), grep( 'KNN', names))
