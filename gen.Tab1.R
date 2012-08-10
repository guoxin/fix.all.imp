beta.set <- exp(seq(log(0.001), log(10), length=30))
lambda.set <- exp(seq(-17, -3))
dyn.load("GK3WN.so")
stringToK3 <- function(string.set, K1, beta){
  len <- length(string.set)
  K3 <- .C("GK3WN", 
    peptides = as.character(string.set), 
    Length = as.integer(len), 
    K3 = as.double(matrix(0.0, nrow = len, ncol = len)), 
    K1 = K1 ^ beta,
    weights = as.double(rep(1, max(nchar(string.set)))), 
    DUP = TRUE, #because of the characters
    PACKAGE = "GK3WN"
   )$K3
   dim(K3) <- c(len, len);gc()
   colnames(K3) <- rownames(K3) <- names(string.set)
   K3
}

getAUC.guoxin <- function(predicted, realValue, threshold){
    AN <- (realValue <= threshold)
    AP <- (!AN)
    compare <- matrix(predicted[AP], ncol = sum(AN), nrow = sum(AP))
    compare <- sweep(x = compare, MARGIN = 2, 
        STATS = predicted[AN], FUN = "-", check.margin = TRUE)
    auc0 <- mean(compare > 0)
    return(max(auc0, 1 - auc0))
}
getRMSE <- function(predicted, realValue){
    return(sqrt(mean( (predicted - realValue)^2 )))
}
getAccuracy.guoxin <- function(predicted, realValue, threshold){
    return( mean((predicted > threshold) == 
                 (realValue > threshold)) )
}

gpgx.time.stamp <- function(...){
    write(print(paste( system("date", intern = TRUE),":",  ...)), 
        file = "performanceRec.txt", append = TRUE)
}

LooReport <- function(pepStrings, K1, rValues, beta.set, lambda.set, threshold, cvInds){
    #Warning: assume "diag" works well below.
    thisReport <- list()
    foldNames <- c("fold1", "fold2", "fold3", "fold4", "fold5")
    rmse.iniMatrix <- matrix(0, ncol = length(lambda.set), nrow = length(beta.set))
    rownames(rmse.iniMatrix) <- beta.set; colnames(rmse.iniMatrix) <- log(lambda.set)
    thisReport$predict <- double(length(pepStrings))
    for(fold.id in 1:5){
        trn.ind <- (cvInds != fold.id); tst.ind <- !trn.ind
        thisReport[[ foldNames[fold.id] ]]$predict <- 
            array(dim = c(sum(tst.ind), length(beta.set), length(lambda.set)))
        thisReport[[ foldNames[fold.id] ]]$rmse.trn <- rmse.iniMatrix
    }
    for(beta.id in 1:length(beta.set)){
        K3Hat <- stringToK3(string.set = pepStrings, K1 = K1, beta = beta.set[beta.id])
        for(fold.id in 1:5){
            trn.ind <- (cvInds != fold.id); tst.ind <- !trn.ind
            decomp <- svd(K3Hat[trn.ind, trn.ind], nv = 0)
            for(lambda.id in 1:length(lambda.set)){
                G.inv <- decomp$u %*% 
    	            sweep(x = t(decomp$u), MARGIN = 1, STATS = 1 / (decomp$d +
    	            sum(trn.ind) * lambda.set[lambda.id]), FUN = "*")
                coeff <- G.inv %*% rValues[trn.ind]
                thisReport[[ foldNames[fold.id] ]]$predict[ ,beta.id, lambda.id] <- 
                    K3Hat[tst.ind, trn.ind] %*% coeff
                thisReport[[ foldNames[fold.id] ]]$rmse.trn[beta.id, lambda.id] <- 
                    sqrt(mean((coeff / diag(G.inv))^2))
            }
        }
    }
    for(fold.id in 1:5){
        bestRMSE <- min(thisReport[[ foldNames[fold.id] ]]$rmse.trn)
        best.ind <- which(thisReport[[ foldNames[fold.id] ]]$rmse.trn == bestRMSE, arr.ind = T)[1, ]
        thisReport[[ foldNames[fold.id] ]]$best.beta.id <- best.ind[1]
        thisReport[[ foldNames[fold.id] ]]$best.lambda.id <- best.ind[2]
        thisReport$predict[cvInds == fold.id] <- 
            thisReport[[ foldNames[fold.id] ]]$predict[ ,best.ind[1], best.ind[2]]
        thisReport[[ foldNames[fold.id] ]]$predict <- NULL
    }
    thisReport$auc <- getAUC.guoxin(predicted = thisReport$predict,
        realValue = rValues, threshold = threshold)
    thisReport$rmse <- getRMSE(predicted = thisReport$predict, realValue = rValues)
    gc(); return(thisReport)
}
load("BLOSUM62.2.Rdata")
load("fixAlleleData")
reports <- list()
for(nm in names(fixAlleleData)){
    gpgx.time.stamp("now computing", nm)
    reports[[nm]] <- LooReport(pepStrings = fixAlleleData[[nm]]$peptides,
        K1 = BLOSUM62.2, rValues = fixAlleleData[[nm]]$rValues,
        beta.set = beta.set, lambda.set = lambda.set, threshold = 0.426, 
        cvInds = fixAlleleData[[nm]]$cvInds)
    save(reports, file = "reports.Rdata")
}
