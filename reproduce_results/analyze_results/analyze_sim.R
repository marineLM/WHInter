source("./analyze_results/compareMethods.R")
# path to the result files obtained by running WHInter, zetaIL, and SPP.
path_to_results <- "./results/"
# path where to save the figures obtained based on these results.
path_to_figures <- "./figures/"
if(!dir.exists(path_to_figures)){
    dir.create(path_to_figures)
}

files_gen <- c("WHInter/Bernoulli_n%d_p%d_qunif_coefnormal_rs0_nnzd_nlam100_rat0.01_max150_Bi1_M0_Bo2_F50_eps1e-08",
               "WHInter/Bernoulli_n%d_p%d_qunif_coefnormal_rs0_nnzd_nlam100_rat0.01_max150_Bi1_M1_Bo2_F50_eps1e-08",
               "WHInter/Bernoulli_n%d_p%d_qunif_coefnormal_rs0_nnzd_nlam100_rat0.01_max150_Bi1_M2_Bo2_F50_eps1e-08",
               "WHInter/Bernoulli_n%d_p%d_qunif_coefnormal_rs0_nnzd_nlam100_rat0.01_max150_Bi1_M2_Bo1_F50_eps1e-08",
               "WHInter/Bernoulli_n%d_p%d_qunif_coefnormal_rs0_nnzd_nlam100_rat0.01_max150_Bi1_M2_Bo0_F50_eps1e-08",
               "zetaIL/Bernoulli_n%d_p%d_qunif_coefnormal_rs0_nnzd_nlam100_rat0.01_max150_Bi1_M0_F50_eps1e-08",
               "zetaIL/Bernoulli_n%d_p%d_qunif_coefnormal_rs0_nnzd_nlam100_rat0.01_max150_Bi1_M1_F50_eps1e-08",
               "zetaIL/Bernoulli_n%d_p%d_qunif_coefnormal_rs0_nnzd_nlam100_rat0.01_max150_Bi1_M2_F50_eps1e-08",
               "SPP/Bernoulli_n%d_p%d_qunif_coefnormal_rs0_nnzd_nlam100_rat0.01_max150_Bi1_F50_eps1e-08")

leg <- c(expression(paste(eta, "(", alpha[l2], ")", "+naive", sep="")),
         expression(paste(eta, "(", alpha[l2], ")", "+MIPS", sep="")),
         expression(paste(eta, "(", alpha[l2], ")", "+IL", sep="")),
         expression(paste(eta, "(", alpha^"*", ")", "+IL", sep="")),
         expression(paste(eta, "(", alpha, "=1", ")", "+IL", sep="")),
         expression(paste(zeta, "+naive", sep="")),
         expression(paste(zeta, "+MIPS", sep="")),
         expression(paste(zeta, "+IL", sep="")),
         "SPP")

leg_bound <- c(expression(paste(eta, "(", alpha[l2], ")", sep="")),
               expression(paste(eta, "(", alpha^"*", ")", sep="")),
               expression(paste(eta, "(", alpha, "=1", ")", sep="")),
               expression(paste(zeta, sep="")),
               "SPP")

col <- c("red", "red", "red", "orange", "green", "deepskyblue", "deepskyblue",  "deepskyblue", "purple")
col_bound <- c("red", "orange", "green", "deepskyblue", "purple")

lty <- c(3, 2, 1, 1, 1, 3, 2, 3, 3)
ind = c(3, 4, 5, 8, 9, 10)
ind_noBlitz = c(3, 4, 5, 8, 9)
leg_short <- c(expression(paste("WHINter - ", eta[alpha[l2]], sep="")),
               expression(paste("WHINter - ", eta[min], sep="")),
               expression(paste("WHINter - ", eta[1], sep="")),
               expression(paste(zeta, "+IL", sep="")),
               "SPP",
               "BLITZ")
leg_noBlitz <- c(expression(paste("WHINter - ", eta[alpha[l2]], sep="")),
               expression(paste("WHINter - ", eta[min], sep="")),
               expression(paste("WHINter - ", eta[1], sep="")),
               expression(paste(zeta, "+IL", sep="")),
               "SPP")


# n fixed p varied
###################
n <- 1000
nfeatures <- c(1000, 3000, 10000)

# Construct the matrix of file names.
files <- matrix("", nrow=length(files_gen), ncol=length(nfeatures))
c = 0
for (p in nfeatures){
    c = c+1
    files[,c] <- sapply(files_gen, function(file){sprintf(file, n, p)})
}
cm1 <- compareMethods2(path_to_results, files)

# Plot times
pdf(paste(path_to_figures, 'WHInter_p.pdf', sep=""), height=4.9, width=4.9)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
matplot(nfeatures, t(cm1$total_time + cm1$preproc_time)/1000, pch=1, ylab = "time (s) (log scale)", xlab = "p (log scale)", log="xy", col = col, type='b', lty=lty, lwd = 1.5)
grid()
legend("topleft", legend=leg, pch=1, col = col, lty=lty, cex=0.8, bty = "n")
dev.off()

# mat_diff records, for each lambda, the sum of squares difference between the vector of coeficients
# obtained with a method against the ones obtained with all other methods.
for (mat_diff in cm1$l_diff){
    print(max(mat_diff))
}

# Bar plots to compare branch check, MIPS, and coordinate descent times
# tmp = t(data.frame(cm1$check_branch_time[, 3], cm1$MIPS_time[, 3], cm1$CD_time[1:8, 3]))
# barplot(as.matrix(tmp))

# p fixed n varied
###################
nsamples <- c(300, 1000, 10000)
p <- 1000
files <- matrix("", nrow=length(files_gen), ncol=length(nsamples))
c = 0
for (n in nsamples){
    c = c+1
    files[,c] <- sapply(files_gen, function(file){sprintf(file, n, p)})
}
cm2 <- compareMethods2(path_to_results, files)

# Plot times
pdf(paste(path_to_figures, 'WHInter_n.pdf', sep=""), height=4.9, width=4.9)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
matplot(nsamples, t(cm2$total_time + cm2$preproc_time)/1000, pch=1, ylab = "time (s) (log scale)", xlab = "n (log scale)", log="xy", col = col, type='b', lty=lty, lwd = 1.5)
grid()
legend("topleft", legend=leg, pch=1, col = col, lty=lty, cex=0.8, bty = "n")
dev.off()

# mat_diff records, for each lambda, the sum of squares difference between the vector of coeficients
# obtained with a method against the ones obtained with all other methods.
for (mat_diff in cm2$l_diff){
    print(max(mat_diff))
}
# Bar plots to compare branch check, MIPS, and coordinate descent times
# tmp = t(data.frame(cm2$check_branch_time[, 1], cm2$MIPS_time[, 1], cm2$CD_time[1:8, 1]))
# barplot(as.matrix(tmp))


# Plot the number of branches opened
####################################
pdf(paste(path_to_figures, '/branches_opened_n300_p1000.pdf', sep=""), height=3.9, width=3.9)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
matplot(data.frame(cm2$branch_opened[,1])[,ind_noBlitz], pch=1, col=col_bound, xlab="lambda index",
        ylab="# violating branches at the first iteration")
grid()
legend(0, 460, legend=leg_noBlitz, pch=1, col = col_bound, cex=0.8, bty = "n")
dev.off()

pdf(paste(path_to_figures, '/branches_opened_n10000_p1000.pdf', sep=""), height=3.9, width=3.9)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
matplot(data.frame(cm2$branch_opened[,3])[,ind_noBlitz], pch=1, col=col_bound, xlab="lambda index",
        ylab="# violating branches at the first iteration")
grid()
legend(-2, 580, legend=leg_noBlitz, pch=1, col = col_bound, cex=0.8, bty = "n")
dev.off()

pdf(paste(path_to_figures, '/branches_opened_n1000_p1000.pdf', sep=""), height=3.9, width=3.9)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
matplot(data.frame(cm2$branch_opened[,2])[,ind_noBlitz], pch=1, col=col_bound, xlab="lambda index",
        ylab="# violating branches at the first iteration")
grid()
legend(0, 550, legend=leg_noBlitz, pch=1, col = col_bound, cex=0.8, bty = "n")
dev.off()


# Sanity checks to see whether the exact same number of branches have been opened for the same bounds
for (i in seq(3)){
    print(all.equal(as.matrix(cm1$branch_opened[,i])[1,1], as.matrix(cm1$branch_opened[,i])[2,1]))
    print(all.equal(as.matrix(cm1$branch_opened[,i])[1,1], as.matrix(cm1$branch_opened[,i])[3,1]))
    print(all.equal(as.matrix(cm1$branch_opened[,i])[6,1], as.matrix(cm1$branch_opened[,i])[7,1]))
    print(all.equal(as.matrix(cm1$branch_opened[,i])[6,1], as.matrix(cm1$branch_opened[,i])[8,1]))
}


# Also compare to Blitz
#######################
p <- 1000
total_time_Blitz_n <- c()
preproc_time_Blitz_n <- c()
for (n in c(300, 1000, 10000)){
    print(n)
    tmp_stats <- read.csv2(paste(path_to_results, "/Blitz/Bernoulli_n", n, "_p", p, "_qunif_coefnormal_rs0_nnzd_nlam100_rat0.01_max150_Bi1_tol1e-08", "_stats.csv", sep=""),
                           header=T, sep=",", dec = ".")
    tmp_preproc <- read.csv2(paste(path_to_results, "/Blitz/Bernoulli_n", n, "_p", p, "_qunif_coefnormal_rs0_nnzd_nlam100_rat0.01_max150_Bi1_tol1e-08", "_preproc.csv", sep=""),
                             header=T, sep=",", dec = ".")
    total_time_Blitz_n <- c(total_time_Blitz_n, sum(tmp_stats$time))
    preproc_time_Blitz_n <- c(preproc_time_Blitz_n, tmp_preproc[[1]])
    
    # See if the model learned by BLITZ is identical to
    # the model learned by WHInter
    tmp_model_B <- read.csv2(paste(path_to_results, "/Blitz/Bernoulli_n", n, "_p", p, "_qunif_coefnormal_rs0_nnzd_nlam100_rat0.01_max150_Bi1_tol1e-08", "_model.csv", sep=""),
                             header=F, sep=",", dec = ".")
    tmp_model_W <- read.csv2(paste(path_to_results, "/WHInter/Bernoulli_n", n, "_p", p, "_qunif_coefnormal_rs0_nnzd_nlam100_rat0.01_max150_Bi1_M2_Bo2_F50_eps1e-08", "_model.csv", sep=""),
                             header=F, sep=",", dec = ".")
    
    print(all.equal(tmp_model_B[,i], tmp_model_W[,i]))
}

total_times_Blitz_n = total_time_Blitz_n + preproc_time_Blitz_n


n <- 1000
total_time_Blitz_p <- c()
preproc_time_Blitz_p <- c()
for (p in c(1000, 3000)){
    print(p)
    tmp_stats <- read.csv2(paste(path_to_results, "/Blitz/Bernoulli_n", n, "_p", p, "_qunif_coefnormal_rs0_nnzd_nlam100_rat0.01_max150_Bi1_tol1e-08", "_stats.csv", sep=""),
                           header=T, sep=",", dec = ".")
    tmp_preproc <- read.csv2(paste(path_to_results, "/Blitz/Bernoulli_n", n, "_p", p, "_qunif_coefnormal_rs0_nnzd_nlam100_rat0.01_max150_Bi1_tol1e-08", "_preproc.csv", sep=""),
                             header=T, sep=",", dec = ".")
    total_time_Blitz_p <- c(total_time_Blitz_p, sum(tmp_stats$time))
    preproc_time_Blitz_p <- c(preproc_time_Blitz_p, tmp_preproc[[1]])
    
    # Check if the model learned by BLITZ is identical to
    # the model learned by WHInter
    tmp_model_B <- read.csv2(paste(path_to_results, "/Blitz/Bernoulli_n", n, "_p", p, "_qunif_coefnormal_rs0_nnzd_nlam100_rat0.01_max150_Bi1_tol1e-08", "_model.csv", sep=""),
                             header=F, sep=",", dec = ".")
    tmp_model_W <- read.csv2(paste(path_to_results, "/WHInter/Bernoulli_n", n, "_p", p, "_qunif_coefnormal_rs0_nnzd_nlam100_rat0.01_max150_Bi1_M2_Bo2_F50_eps1e-08", "_model.csv", sep=""),
                             header=F, sep=",", dec = ".")
    
    print(all.equal(tmp_model_B[,i], tmp_model_W[,i]))
}

total_times_Blitz_p = total_time_Blitz_p + preproc_time_Blitz_p



to_plot_p = rbind((cm1$total_time + cm1$preproc_time)/1000, c(total_times_Blitz_p, NA))
pdf(paste(path_to_figures, '/WHInter_p_withBlitz.pdf', sep=""), height=3.9, width=3.9)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
matplot(nfeatures, t(to_plot_p[ind,]), pch=1, ylab = "time (s) (log scale)", xlab = "p (log scale)", log="xy",
        col = c(col, "black")[ind], type='b', lty=c(lty, 3)[ind], lwd = 1.5)
grid()
legend("topleft", legend=leg_short, pch=1, col = c(col, "black")[ind], lty=c(lty, 3)[ind], cex=0.8, bty = "n")
dev.off()


to_plot_n = rbind((cm2$total_time + cm2$preproc_time)/1000, total_times_Blitz_n)
pdf(paste(path_to_figures, '/WHInter_n_withBlitz.pdf', sep=""), height=3.9, width=3.9)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
matplot(nsamples, t(to_plot_n[ind, ]), pch=1, ylab = "time (s) (log scale)", xlab = "n (log scale)", log="xy",
        col = c(col, "black")[ind], type='b', lty=c(lty, 3)[ind], lwd = 1.5)
grid()
legend("topleft", legend=leg_short, pch=1, col = c(col, "black")[ind], lty=c(lty, 3)[ind], cex=0.8, bty = "n")
dev.off()
