## A script executing comparisons between different local models and different prewarpings. 
## Used to create Figs 4-6.

library(laGP)
library(lhs)
library(activegp)
library(hetGP)
library(FNN)

#sourceCpp("~/asl/code/examples/kernelexps.cpp")
source('R/test_functions.R') # Bring in test functions.
source('R/sample_C.R') # The empirical measure estimator.
source('R/cli.R') # The empirical measure estimator.

## Handle input arguments
args <- commandArgs(trailingOnly = TRUE)
ret <- process_cli(args)
PROBLEMS <- ret$PROBLEMS
DR_METHODS <- ret$DR_METHODS

scorep <- function(YY, mu, s2) { mean(-(mu - YY)^2/s2 - log(s2)) }

#N_FUNC <- 40000
N_FUNC <- 400
NN_FUNC <- 2000
prop_test <- 0.1
MIN_DR <- 5
MAX_DR <- 25
BY_DR <- 1
VAL_PROP <- 0.2 
CV_GRID_SIZE <- 8
KNN_K <- 5

#gp_subset_size <- 1500 # Size of the subset GP comparator.
#n_Cs <- 5# Number of subsets to estimate C
#max_subsize <- 1500 # Size of the subsets to estimate C
gp_subset_size <- 300 # Size of the subset GP comparator.
n_Cs <- 5# Number of subsets to estimate C
max_subsize <- 300 # Size of the subsets to estimate C

subbag <- TRUE

scorep <- function(YY, mu, s2) { mean(-(mu - YY)^2/s2 - log(s2)) }

if (subbag) {
    library(doFuture)
    library(foreach)
    library(doRNG)
    registerDoFuture()
    plan(multisession)
}

C_GP_subbag <- function(X, y, n_Cs, subsize, measure = 'lebesgue', theta = NULL) {

    N <- nrow(X)
    subsamps <- lapply(1:n_Cs, function(i) sample(N, subsize))

    fits <- list()
    for (subset in subsamps) {
        if (is.null(theta)) {
            bounds <- hetGP:::auto_bounds(X = X[subset,], covtype = 'Gaussian')
            fits[[length(fits)+1]] <- mleHomGP(X[subset,,drop=F], y[subset], lower = bounds$lower, upper = bounds$upper)
        } else {
            fits[[length(fits)+1]] <- mleHomGP(X[subset,,drop=F], y[subset], known = list(theta))
        }
    }

    Cs <- foreach(fit = fits, .export = c("X", "y", "measure", "theta")) %dorng% {
        library(hetGP)
        library(activegp)
        source('R/sample_C.R') # Bring in test functions.
        if (measure == 'empirical')
            {return(C_GP_empirical(fit))}
        else if (measure == 'lebesgue')
            {return(C_GP(fit)$mat)}
        else {
            stop("Unsupported Measure.")
        }
    }

    Ch <- Reduce(function(x,y) x+y, Cs) / n_Cs
    return(Ch)
}


iters <- 1

m1 <- c("aGP", "aGP_trunc", "Vecchia", "Vecchia_trunc", "knn", "knn_trunc")
proj_mnames <- as.character(sapply(DR_METHODS, function(drm) paste(drm, m1, sep = '-')))
method_names <- c("aGP", "subsetGP", "vecchia", "knn", proj_mnames)
ncomp <- length(method_names)
ndr <- length(DR_METHODS)

# Start with test functions
test_probs <- lapply(testfuns, function(x) {x$isfun <- TRUE; x})
test_probs <- Filter(function(prob) prob$name %in% PROBLEMS, test_probs)

# Bring in datasets
if ('CNC' %in% PROBLEMS) {
    cnc <- read.csv('./data/communities.data', header = FALSE)
    cnc <- as.data.frame(apply(cnc, 2, function(x) as.numeric(as.character(x))))
    fileName <- "./data/communities.colnames"
    colnames(cnc) <- strsplit(readChar(fileName, file.info(fileName)$size), split = ', ')[[1]]

    cnc <- cnc[,complete.cases(t(cnc))]
    cnc_X <- cnc[3:(ncol(cnc)-1)]
    cnc_y <- cnc$ViolentCrimesPerPop
    test_probs$cnc <- list(X = as.matrix(cnc_X), y = cnc_y, isfun = FALSE, name = "CNC")
}
if ('Temperature' %in% PROBLEMS) {
    tX <- read.table('./data/temp_train.inputs', header = FALSE)
    tX <- (tX - min(tX)) / (max(tX) - min(tX))
    ty <- as.numeric(as.matrix(read.table('./data/temp_train.targets', header = FALSE)))
    test_probs$temp <- list(X = as.matrix(tX), y = ty, isfun = FALSE, name = "Temperature")
}
if ('CT Scans' %in% PROBLEMS) {
    ct_dat <- read.csv('./data/slice_localization_data.csv')
    tX <- ct_dat[,grep('value', colnames(ct_dat))]
    tX <- (tX - min(tX)) / (max(tX) - min(tX))
    ty <- ct_dat$reference
    test_probs$ct <- list(X = as.matrix(tX), y = ty, isfun = FALSE, name = "CT Scans")
}
if ('MOPTA' %in% PROBLEMS) {
    #mopta <- read.csv('./data/random_mopta.csv')
    mopta <- read.csv('/storage/nate/mopta_500k.csv')
    tX <- mopta[2:125]
    #tX <- (tX - min(tX)) / (max(tX) - min(tX)) # Already in [0,1]
    Ys <- mopta[,126:ncol(mopta)]
    Ys <- apply(Ys, 2, function(y) (y - mean(y)) / sd(y))
    ty <- rowMeans(Ys)
    test_probs$mopta <- list(X = as.matrix(tX), y = ty, isfun = FALSE, name = "MOPTA")
}

NF <- length(test_probs)

tt <- Sys.time()

mses <- array(NA, dim = c(NF, iters, ncomp))
scores <- array(NA, dim = c(NF, iters, ncomp))
times <- array(NA, dim = c(NF, iters, ncomp))
c_times <- array(NA, dim = c(NF, iters, ndr))
dr_outs <- array(NA, dim = c(NF, iters, ndr)) # Stores truncated dimension
for (f_i in 1:length(test_probs)) {
    fobj <- test_probs[[f_i]]
    if (fobj$isfun) {
        target_func <- fobj$fun
        m <- fobj$d
        n <- N_FUNC
        nn <- NN_FUNC
    } else {
        m <- ncol(fobj$X)
        n <- round(nrow(fobj$X)*(1-prop_test))
        nn <- nrow(fobj$X) - n
    }
    cur_max_dr <- pmin(MAX_DR, m)

    subsize <- pmin(n, max_subsize)

    for (iter in 1:iters) {
        print(iter)
        if (fobj$isfun) {
            x <- randomLHS(n + nn, m)
            y <- apply(x, 1, target_func)
            X <- x[1:n,]
            Y <- y[1:n]
            XX <- x[-(1:n),]
            YY <- y[-(1:n)]
        } else {
            train_inds <- sample(nrow(fobj$X), n)
            X <- fobj$X[train_inds,]
            Y <- fobj$y[train_inds]
            XX <- fobj$X[-train_inds,]
            YY <- fobj$y[-train_inds]
        }

        ## Get validation set for choosing truncation size.
        valset <- sample(n, round(VAL_PROP*n))

        # LAgp 
        omp_nth <- 40
        tt_aGP <- Sys.time()

        mod_aGP <- aGPsep(X, Y, XX, omp.threads=omp_nth, verb=0)

        scores[f_i,iter,1] <- scorep(YY, mod_aGP$mean, mod_aGP$var)
        mses[f_i,iter,1] <- mean((YY - mod_aGP$mean)^2)
        times[f_i,iter,1] <- as.numeric(Sys.time() - tt_aGP, unit = 'secs')

        # Subset GP Regular 
        tt_subGP <- Sys.time()
        nsub <- pmin(gp_subset_size, n)
        d2 <- darg(list(mle=TRUE, max=100), X)
        subs <- sample(1:nrow(X), nsub, replace=FALSE)
        gpsi <- newGPsep(X[subs,], Y[subs], rep(d2$start, m), g=1/1000, dK=TRUE)
        that <- mleGPsep(gpsi, tmin=d2$min, tmax=d2$max, ab=d2$ab, maxit=200)
        pred_subGP <- predGPsep(gpsi, XX, lite=TRUE)
        deleteGPsep(gpsi)

        scores[f_i,iter,2] <- scorep(YY, pred_subGP$mean, pred_subGP$s2)
        mses[f_i,iter,2] <- mean((YY - pred_subGP$mean)^2)
        times[f_i,iter,2] <- as.numeric(Sys.time() - tt_subGP, unit = 'secs')

        if (m < 10) {
            tt_vecc <- Sys.time()
            fit.rv <- GpGp::fit_model(Y,X,covfun_name = 'exponential_isotropic', silent = TRUE)
            preds <- GpGp::cond_sim(fit=fit.rv,locs_pred=XX, X_pred=matrix(1,nrow=nrow(XX),ncol=1), nsims = 30)
            pred_mean <- rowMeans(preds)
            pred_var <- apply(preds, 1, var)

            scores[f_i,iter,3] <- scorep(YY, pred_mean, pred_var)
            mses[f_i,iter,3] <- mean((YY - pred_mean)^2)
            times[f_i,iter,3] <- as.numeric(Sys.time() - tt_vecc, unit = 'secs')
        } else {
            scores[f_i,iter,3] <- NA
            mses[f_i,iter,3] <- NA
            times[f_i,iter,3] <- NA
        }

        # KNN
        tt_knn <- Sys.time()
        preds <- knn.reg(X, y = Y, test = XX)$pred

        ind <- 4
        scores[f_i,iter,ind] <- NA
        mses[f_i,iter,ind] <- mean((YY - preds)^2)
        times[f_i,iter,ind] <- as.numeric(Sys.time() - tt_knn, unit = 'secs')


        dr_range <- seq(from = MIN_DR, to = cur_max_dr, by = BY_DR)

        ########################################################################################
        for (dr_method in DR_METHODS) {
            tt_C <- Sys.time()
            dr_ind <- which(DR_METHODS==dr_method)
            if (dr_method == 'lengthscale') {
                # Subset GP Regular 
                tt_subGP <- Sys.time()

                tt_agps <- Sys.time()
                scale <- sqrt(that$d)
                L <- t(diag(1/scale)[rank(scale),])

            } else if (dr_method == 'lebesgue') {
                # Rotated laGP -- Lebesgue
                if (subbag) {
                    Ch <- C_GP_subbag(X, Y, n_Cs, subsize, measure = 'lebesgue')
                } else {
                    psubhetGP <- mleHomGP(X = X[subs,], Y[subs])
                    Ch <- C_GP(psubhetGP)$mat
                }
                Ch <- Ch / sum(diag(Ch))
                Ched <- eigen(Ch)
                L <- Ched$vectors %*% diag(sqrt(Ched$values))
            } else if (dr_method == 'empirical') {
                if (subbag) {
                    Ch <- C_GP_subbag(X, Y, n_Cs, subsize, measure = 'empirical')
                } else {
                    stop("Only subbagging for empirical.")
                }
                Ch <- Ch / sum(diag(Ch))
                Ched <- eigen(Ch)
                L <- Ched$vectors %*% diag(sqrt(Ched$values))
            } else if (dr_method == 'identity') {
                L <- diag(m)
            } else if (substr(dr_method,1,4) == 'iso_') {
                drm <- substring(dr_method, 5)
                subsamp <- sample(n, nsub)
                dists <- as.matrix(dist(X[subsamp,]))
                sigma_med <- median(dists[upper.tri(dists)])
                l_cands <- 2*seq(from = 0.5 * sigma_med, to = CV_GRID_SIZE * sigma_med) #Factor of 2 adjusts for different kernel definitions.

                knn_err <- rep(NA, CV_GRID_SIZE)
                Ls <- list()
                for (theta in l_cands) {
                    Ct <- C_GP_subbag(X, y, n_Cs, subsize, measure = drm, theta = theta)
                    Ct <- Ct / sum(diag(Ct))
                    Cted <- eigen(Ct)
                    L <- Cted$vectors %*% diag(sqrt(Cted$values))
                    Ls[[length(Ls)+1]] <- L
                    Xr <- X %*% L
                    knn_err[which(theta==l_cands)] <- sum(knn.reg(Xr, y = Y)$residuals^2)
                }
                L <- Ls[[which.min(knn_err)]]
            } else {
                stop(paste("Unsupported dr_method: ", dr_method, collapse = ''))
            }

            Xr <- X %*% L
            XXr <- XX %*% L


            # Pick dimension reduction size. 
            knn_err <- rep(NA, length(dr_range))
            for (dr in dr_range) {
                Xt <- Xr[,1:dr,drop=FALSE]
                knn_err[which(dr_range==dr)] <- sum(knn.reg(Xt, y = Y)$residuals^2)
            }
            bic <- length(subs)*log(knn_err/length(subs)) + (MIN_DR:cur_max_dr)*log(length(subs))
            dr_out <- (MIN_DR:cur_max_dr)[which.min(bic)]
            dr_outs[f_i,iter,dr_ind] <- dr_out
            Xt <- Xr[,1:dr_out,drop=FALSE]
            XXt <- XXr[,1:dr_out,drop=FALSE]

            c_times[f_i,iter,dr_ind] <- as.numeric(Sys.time() - tt_C, unit = 'secs')

            ## later: PCA
            tt_agpr <- Sys.time()
            mod_aGPr  <- aGP(Xr, Y, XXr, omp.threads=omp_nth, verb=0)

            ind <- 4+length(m1)*(dr_ind-1) + 1
            scores[f_i,iter,ind] <- scorep(YY, mod_aGPr$mean, mod_aGPr$var)
            mses[f_i,iter,ind] <- mean((YY - mod_aGPr$mean)^2)
            times[f_i,iter,ind] <- as.numeric(Sys.time() - tt_agpr, unit = 'secs')

            ## Truncated laGP
            tt_agpt <- Sys.time()
            mod_aGPt <- aGP(Xt, Y, XXt, omp.threads=omp_nth, verb=0)

            ind <- 4+length(m1)*(dr_ind-1) + 2
            scores[f_i,iter,ind] <- scorep(YY, mod_aGPt$mean, mod_aGPt$var)
            mses[f_i,iter,ind] <- mean((YY - mod_aGPt$mean)^2)
            times[f_i,iter,ind] <- as.numeric(Sys.time() - tt_agpt, unit = 'secs')

            # Vecchia
            ind <- 4+length(m1)*(dr_ind-1) + 3
            if (m < 10) {
                tt_vecc <- Sys.time()
                fit.rv <- GpGp::fit_model(Y,Xr,covfun_name = 'exponential_isotropic', silent = TRUE)
                preds <- GpGp::cond_sim(fit=fit.rv,locs_pred=XXr, X_pred=matrix(1,nrow=nrow(XXr),ncol=1), nsims = 30)
                pred_mean <- rowMeans(preds)
                pred_var <- apply(preds, 1, var)

                scores[f_i,iter,ind] <- scorep(YY, pred_mean, pred_var)
                mses[f_i,iter,ind] <- mean((YY - pred_mean)^2)
                times[f_i,iter,ind] <- as.numeric(Sys.time() - tt_vecc, unit = 'secs')
            } else {
                scores[f_i,iter,ind] <- NA
                mses[f_i,iter,ind] <- NA
                times[f_i,iter,ind] <- NA
            }

            ## Truncated Vecchia
            tt_vgpt <- Sys.time()
            fit.rv=GpGp::fit_model(Y,Xt, silent = TRUE, covfun_name = 'exponential_isotropic')
            preds <- GpGp::cond_sim(fit=fit.rv,locs_pred=XXt, X_pred=matrix(1,nrow=nrow(XXr),ncol=1), nsims = 30)
            pred_mean <- rowMeans(preds)
            pred_var <- apply(preds, 1, var)

            ind <- 4+length(m1)*(dr_ind-1) + 4
            scores[f_i,iter,ind] <- scorep(YY, pred_mean, pred_var)
            mses[f_i,iter,ind] <- mean((YY - pred_mean)^2)
            times[f_i,iter,ind] <- as.numeric(Sys.time() - tt_vgpt, unit = 'secs')

            # KNN
            tt_knn <- Sys.time()
            preds <- knn.reg(Xr, y = Y, test = XXr)$pred

            ind <- 4+length(m1)*(dr_ind-1) + 5
            scores[f_i,iter,ind] <- NA
            mses[f_i,iter,ind] <- mean((YY - preds)^2)
            times[f_i,iter,ind] <- as.numeric(Sys.time() - tt_knn, unit = 'secs')

            ## Truncated KNN
            tt_knnt <- Sys.time()
            preds <- knn.reg(Xt, y = Y, test = XXt)$pred

            ind <- 4+length(m1)*(dr_ind-1) + 6
            scores[f_i,iter,ind] <- NA
            mses[f_i,iter,ind] <- mean((YY - preds)^2)
            times[f_i,iter,ind] <- as.numeric(Sys.time() - tt_knnt, unit = 'secs')

        }
    }
}



test_prob_names <- sapply(test_probs, function(tp) tp$name)
nice_time <- gsub(' ', '_', as.character(Sys.time()))
title <- paste("C_", n_Cs, "_subsize_", max_subsize, "_at_", nice_time, '.RData', sep = '')
save(mses, scores, dr_outs, times, c_times, method_names, test_prob_names, file = file.path('data', 'runs', title))

print("Total Time:")
print(Sys.time() - tt)
print("Sum of recorded times:")
print(sum(times, na.rm = T) + sum(c_times))
