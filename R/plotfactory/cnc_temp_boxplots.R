#!/usr/bin/Rscript
#  R/make_cnc_temp_boxplots.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 12.08.2020

load("./data/runs/latest_ed.RData")
source("R/plotfactory/common.R")

scalemse <- length(names)/length(hasscore)*0.9# 0.9 is Empirical adjustment for margins based on on-screen differences.

boldnames <- as.expression(rep(NA,length(names)))
tobold <- grep('^.-', names)
for (i in 1:length(names)) {
    if (i %in% tobold) {
        #names[i] <- as.expression(paste('bold(', names[[i]], ')', sep = ''))
        boldnames[i] <- as.expression(bquote(bold(.(names[[i]]))))
    } else {
        boldnames[i] <- as.expression(bquote(.(names[[i]])))
    }
}

pdf("images/cnc.pdf", width = 8, height = 4*height_scale)
par(oma = oma, mar = mar)
layout(matrix(1:2, nrow = 1), widths = c(scalemse,1))
title <-  "Communities and Crime"
f_i <- 1
boxplot(log10(mses[f_i,,]), names = boldnames, las = 2, ylab = "logMSE", col = cols, border = slightly_darker_cols, cex.lab = cex_axis)
boxplot(-scores[f_i,,hasscore], names = boldnames[hasscore], las = 2, ylab = "-Score", col = cols[hasscore], border = slightly_darker_cols[hasscore], cex.lab = cex_axis)
mtext(title, outer = TRUE, cex = 1.5*height_scale, adj = 1)
dev.off()

pdf("images/temperature.pdf", width = 8, height = 4*height_scale)
par(oma = oma, mar = mar)
layout(matrix(1:2, nrow = 1), widths = c(scalemse,1))
title <-  "Temperature"
f_i <- 2
boxplot(log10(mses[f_i,,]), names = boldnames, las = 2, ylab = "logMSE", col = cols, border = slightly_darker_cols, cex.lab = cex_axis)
boxplot(-scores[f_i,,hasscore], names = boldnames[hasscore], las = 2, ylab = "-Score", col = cols[hasscore], border = slightly_darker_cols[hasscore], cex.lab = cex_axis)
mtext(title, outer = TRUE, cex = 1.5*height_scale, adj = 1)
dev.off()
