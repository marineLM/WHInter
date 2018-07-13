compareMethods2 <- function(path_to_results, files){
    
    nmeth = nrow(files) # Number of methods compared
    ndata = ncol(files) # Number of datasets studied
    
    total_time <- data.frame(matrix(nrow=nmeth, ncol=ndata))
    preproc_time <- data.frame(matrix(nrow=nmeth, ncol=ndata))
    check_branch_time <- data.frame(matrix(nrow=nmeth, ncol=ndata))
    MIPS_time <- data.frame(matrix(nrow=nmeth, ncol=ndata))
    CD_time <- data.frame(matrix(nrow=nmeth, ncol=ndata))
    proj_time <- data.frame(matrix(nrow=nmeth, ncol=ndata))
    safeprune1_time <- data.frame(matrix(nrow=nmeth, ncol=ndata))
    branch_opened <- matrix(list(), nrow=nmeth, ncol=ndata)
    activeset <- matrix(list(), nrow=nmeth, ncol=ndata)
    features_selected <- matrix(list(), nrow=nmeth, ncol=ndata)
    l_diff <- list()
    
    for (j in seq(ndata)){
        l_model <- list()
        for (i in seq(nmeth)){
            tmp_stats <- read.csv2(paste(path_to_results, files[i, j], "_stats.csv", sep=""), header=T, sep=" ", dec = ".")
            total_time[i, j] <- sum(tmp_stats$time)
            features_selected[i, j] <- list(tmp_stats$n_features_selected)
            CD_time[i, j] <- sum(tmp_stats[, grepl("time_CD", colnames(tmp_stats))])
            if ( substring(files[i, j], 1, 3) != "SPP"){
                branch_opened[i, j] <- list(tmp_stats$n_branch_opened1)
                activeset[i, j] <- list(tmp_stats$n_activeset1)
                check_branch_time[i, j] <- sum(tmp_stats[, grepl("time_check", colnames(tmp_stats))])
                MIPS_time[i, j] <- sum(tmp_stats[, grepl("time_MIPS", colnames(tmp_stats))])
            }else{
                branch_opened[i, j] <- list(tmp_stats$n_branch_opened_safepruning)
                activeset[i, j] <- list(tmp_stats$n_safeset)
                proj_time[i, j] <- sum(tmp_stats$time_proj1)
                safeprune1_time[i, j] <- sum(tmp_stats$time_safeprune)
            }
            tmp_preproc <- read.csv2(paste(path_to_results, files[i, j], "_preproc.csv", sep=""), header=T, sep=",", dec = ".")
            preproc_time[i, j] <- tmp_preproc$time_preproc

            # Check if learned models are identical
            tmp_model <- read.csv2(paste(path_to_results, files[i, j], "_model.csv", sep=""), header=F, sep=",", dec = ".")
            l_model <- c(l_model, list(tmp_model))
        }

        ilammax = l_model[[1]][nrow(l_model[[1]]),1]
        mat_diff <- matrix(0, ilammax, nmeth)
        for (ilam in seq(ilammax)){
            print(ilam)
            l_model_ilam <- lapply(l_model, function(tmp_model_i){
                ind <- which(tmp_model_i[,1] == ilam)
                tmp_model_sub <- tmp_model_i[ind, ]
                ord <- order(tmp_model_sub[, 5])
                tmp_model_sub <- tmp_model_sub[ord, ]
                is_greater <- tmp_model_sub[, 3] > tmp_model_sub[, 4]
                tmp <- tmp_model_sub[is_greater, 3]
                tmp_model_sub[is_greater, 3] <- tmp_model_sub[is_greater, 4]
                tmp_model_sub[is_greater, 4] <- tmp
                return(tmp_model_sub)
            })
            
            v_length <- sapply(l_model_ilam, nrow)
            if (max(v_length)!=min(v_length)){
                message("The number of selected features differs according to the method used for the ",
                        ilam, "th lambda and ndata/nmeth=", j, "/", i)
                print(v_length)
            }else{
                # Check if the features match
                ref_model_ilam <- l_model_ilam[[1]]
                l_differences <- lapply(l_model_ilam, function(tmp_model_ilam_i){
                    which(tmp_model_ilam_i[ ,3:4] != ref_model_ilam[, 3:4], arr.ind=T)
                })
                for (i in seq(length(l_differences))){
                    if(nrow(l_differences[[i]]) != 0){
                        diff_ind <- l_differences[[i]][, 1]
                        message("Selected features IDs differ between method 1 and the ", i, "th method.")
                        message("The feature(s) for which there is a problem is (are) associated to coef ")
                        print(l_model_ilam[[i]][diff_ind , 3:5])
                        message("and")
                        print(ref_model_ilam[diff_ind , 3:5])
                    }
                }
                
                # Check if the coef match
                coef_diff <- sapply(l_model_ilam, function(tmp_model_ilam_i){
                    s = sum((ref_model_ilam[,5] - tmp_model_ilam_i[, 5])**2)
                    return(s)
                })
                
                mat_diff[ilam,] <- coef_diff
            }
        }
        l_diff <- append(l_diff, list(mat_diff))
    }
    return(list(total_time=total_time, preproc_time=preproc_time, check_branch_time=check_branch_time,
                MIPS_time=MIPS_time, CD_time=CD_time, proj_time=proj_time,
                safeprune1_time=safeprune1_time, branch_opened=branch_opened,
                activeset=activeset, features_selected=features_selected, l_diff=l_diff))
}