path_to_data <- "./data_preprocessed/"
path_to_results <- "./results/"
path_to_figures <- "./figures/"
if(!dir.exists(path_to_figures)){
    dir.create(path_to_figures)
}

l_recovery <- list()

dims <- data.frame(c(1000, 1000), c(1000, 3000), c(1000, 10000), c(300, 1000), c(10000, 1000))

pdf(paste(path_to_figures, 'support_recovery.pdf', sep=""), height=3.5, width=5)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
plot.new()
# The indexes and coefficients of the features in the support are the same for all combinations of n and p
gt <- read.csv2(paste(path_to_data, "GroundTruth_n", 1000, "_p", 1000, "_qunif_coefnormal_rs0_nnzd.csv", sep=""), header=F, sep=",", dec = ".")
start = round(min(gt[,3]), 1)
end = round(max(gt[,3]), 1)
by = (end-start)/10
plot.window(xlim=c(start, end), ylim=c(0, 2*length(dims)-1.5))

for (i in seq(ncol(dims))){
    n <- dims[1, i]
    p <- dims[2, i]
    # Get ground truth
    gt <- read.csv2(paste(path_to_data, "GroundTruth_n", n, "_p", p, "_qunif_coefnormal_rs0_nnzd.csv", sep=""), header=F, sep=",", dec = ".")
    for (k in seq(nrow(gt))){
        if(gt[k, 1] > gt[k, 2]){
            tmp <- gt[k, 1]
            gt[k, 1] <- gt[k, 2]
            gt[k, 2] <- tmp
        }
    }
    # Get selected features
    pred <- read.csv2(paste(path_to_results, "WHInter/Bernoulli_n", n, "_p", p, "_qunif_coefnormal_rs0_nnzd_nlam100_rat0.01_max150_Bi1_M2_Bo2_F50_eps1e-08_model.csv", sep=""), header=F, sep=",", dec = ".")
    ilams <- unique(pred[,1])
    recovery <- data.frame(matrix(0, nrow=length(ilams), ncol=4))
    colnames(recovery) <- c('TP', 'total', 'mean_coef', 'sd_coef')
    for (ilam in ilams){
        ind <- which(pred[,1] == ilam)
        pred_ilam <- pred[ind,]
        matches <- apply(pred_ilam[, 3:4], 1, function(row){
            candidate <- min(row) + max(row)*1i
            match <- which(gt[,1] + gt[,2]*1i == candidate)
            return(match)
        })
        matches <- as.numeric(Filter(length, matches))
        recovery$TP[ilam] = length(matches)
        recovery$total[ilam] = nrow(pred_ilam)
        recovery$mean_coef[ilam] = mean(abs(gt[matches, 3]))
        recovery$sd_coef[ilam] =  sd(abs(gt[matches, 3]))
    }
    l_recovery <- append(l_recovery, list(recovery))
    recovered <- rep("black", nrow(gt))
    recovered[matches] <- "red"
    points(gt[,3], rep(2*i-1.5, nrow(gt)), col=recovered)
}

axis(1, at=seq(start, end, by=by), labels=seq(start, end, by=by))
ylab <- c("n=1000\n p=1000", "n=1000\n p=3000", "n=1000\n p=10000",
          "n=300\n p=1000", "n=10000\n p=1000")
axis(2, at=2*seq(length(dims))-1.5, labels=ylab, las=1, tick=F, cex.axis=0.9)
title(xlab="true coefficients of the support")
# legend(-2.7, 2*length(dims)+2, legend=c("selected", "not seleted"), col=c("red", "black"), pch=1,
#        bg="grey", bty="o", box.col="white")
dev.off()

