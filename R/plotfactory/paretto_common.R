

iters <- dim(mses)[2]
names <- rep(names, rep(iters, length(names)))
cols <- rep(cols, rep(iters, length(cols)))
slightly_darker_cols <- rep(slightly_darker_cols, rep(iters, length(slightly_darker_cols)))
locmod <- gsub('-T$', '', gsub('^.-', '', names))

drt <- gsub('-.*', '', names)
drt[nchar(drt) > 1] <- 'N'

usymbols <- list()
udrt <- unique(drt)
for (dr in udrt) {
    usymbols[[dr]] <- (21:25)[which(udrt==dr)]
}
symbols <- unname(unlist(usymbols[drt]))
