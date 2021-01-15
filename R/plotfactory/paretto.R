#!/usr/bin/Rscript
#  /home/nate/rot_vec/R/plotfactory/paretto.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.13.2021

dataf <- c("latest_abc", "latest_ed", "latest_g")
pointalpha <- 0.15

titledict <- list('CNC' = 'Communities and Crime', 'Temperature' = 'Temperature', 'piston' = 'Piston', 'borehole' = 'Borehole', 'MOPTA' = 'MOPTA', 'robot arm' = 'Robot Arm')

for (df in dataf) {
    fname <- file.path("data/runs", paste(df, '.RData', sep = ''))
    load(fname)
    source("R/plotfactory/common.R")
    medmse <- apply(mses, c(1,3), median)
    medscores <- apply(scores, c(1,3), median)
    sel_every <- dim(mses)[2]
    source("R/plotfactory/paretto_common.R")

    for (f_i in 1:dim(mses)[1]) {

        mses_i <- c(mses[f_i,,])
        scores_i <- c(scores[f_i,,])

        medcols <- cols[(1:length(cols)-1)%%10==0]
        #medsdc <- slightly_darker_cols[(1:length(cols)-1)%%10==0]
        medsdc <- 'black'
        slightly_darker_cols <- 'black'
        medsymbs <- symbols[(1:length(cols)-1)%%10==0]

        hastrunc <- grep('-T', names)
        cexvec <- rep(3, length(names))
        cexvec[hastrunc] <- 0
        medcexvec <- cexvec[(1:length(cols)-1)%%10==0]

        savename <- paste('images/paretto_', gsub('\\s','_', test_prob_names[f_i]), '.pdf', sep='')
        pdf(savename, width = 4, height = 4*height_scale)
        par(oma = oma, mar = mar)
        #pdf(savename, width = 8, height = 8)
        c1 <- rgb(t(col2rgb(cols)/255), alpha = pointalpha)
        c2 <- rgb(t(col2rgb(slightly_darker_cols)/255), alpha = pointalpha)
        plot(log10(mses_i), -scores_i, col = c2, bg = c1, pch = symbols, xlab = 'logMSE', ylab = '-Score', cex = 1.5, lwd = cexvec, cex.lab = cex_axis)
        points(log10(medmse[f_i,]), -medscores[f_i,], col = medsdc, bg = medcols, pch = medsymbs, lwd = medcexvec, cex = 1.5)
        #legend('topleft', legend = unique(locmod), fill = unique(cols), col = unique(slightly_darker_cols))
        legend('topleft', legend = gsub('N','',unique(drt)), pch = unique(symbols))
        title <- titledict[[test_prob_names[f_i]]]
        #mtext(title, outer = TRUE, cex = 1.5)
        dev.off()
    }
}

