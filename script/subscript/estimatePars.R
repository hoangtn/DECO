estimatePars <- function(pars, mcmcResult, nThin = NULL, eSD = 10^-8){
    mcmcDataFrame <- as.data.frame(mcmcResult)
        if (!is.null(nThin))
        mcmcDataFrame <- mcmcDataFrame[seq(1, dim(mcmcDataFrame)[1], by = nThin),]
    allPars <- colnames(mcmcDataFrame)
#calculate standard deviations for alpha's 

    
    if (is.null(pars)){
        message("\nYou are estimating all parameters, you can input a specific parameter(s) for pars\n")
        pars = allPars
    }

#    pars <- pars[grep("hyper|pi|alpha", pars)]
    message("====\nOnly pi, alpha and hyper parameters are estimated in this step\n",
            "gTADA does not calculate HPDs for hyper betas, just their medians\n===\n")

    if (length(pars[!is.element(pars, colnames(mcmcDataFrame))]) > 0)
        warning((pars[!is.element(pars, colnames(mcmcDataFrame))]), " is/are not in mcmc chains")
    pars <- pars[is.element(pars, colnames(mcmcDataFrame))]
#print(pars)
    hpdList <- NULL
    bSD <- apply(mcmcDataFrame[, pars], 2, sd)
    
    for (iPar in pars) {
#        message("Estimating for ", iPar)
        xData <-  mcmcDataFrame[, iPar]
        if ((sd(xData) > eSD) & (length(grep("Beta", iPar)) == 0))
            hpdList <- rbind(hpdList, loc1stats(xData))
        else
            hpdList <- rbind(hpdList, rep(median(xData), 3))
    }
    rownames(hpdList) <- pars


    hpdList <- cbind(hpdList, bSD, hpdList[, 1]/bSD)
    hpdList <- cbind(hpdList, pnorm(abs(hpdList[, 1]/bSD), lower.tail = FALSE))
                     #rep(dim(mcmcDataFrame)[1], dim(hpdList)[1]))
    hpdList <- cbind(hpdList)

    colnames(hpdList) <- c("Mode", "lCI", "uCI", "pSD", "zScore", "pValue")

return(hpdList)
}
