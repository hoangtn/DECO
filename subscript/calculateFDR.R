calculateFDR <- function(pars, geneSet = NULL,
                         caseData = NULL,
                         controlData = NULL,
                         dnData = NULL,
                         mutData = NULL,
                         geneName){
    outData <- data.frame(geneName)
    if (!is.null(dnData))
        outData <- cbind(outData, dnData)
    if (!is.null(mutData))
        outData <- cbind(outData, mutData)
    if (!is.null(caseData))
        outData <- cbind(outData, caseData)
    if (!is.null(controlData))
        outData <- cbind(outData, controlData)


    pi0 = pars$alpha0[1]
    if (!is.null(geneSet)){
        geneSet <- data.frame(geneSet)
        for (i in 1:dim(geneSet)[2])
            pi0 = pi0 + pars$alpha0[i+1]*geneSet[, i]
        } else {
            message("\nNo gene set in this calculation\n")
        }
    
    pi0 <- exp(pi0)/(1 + exp(pi0))



    bfAll <- llk1All <- llk0All <- rep(1, dim(outData)[1])

    if ( length(pars$gammaMeanDN) == 0) {
        message("No parameters for de novo data; therefore, these categories are not calculated in this step.\n")
        }  else {
        bfDN <- llk1DN <- llk0DN <- matrix(1, nrow = dim(dnData)[1], ncol = dim(dnData)[2]) ##De novo bayes factors
        for (j2 in 1:dim(bfDN)[2]) {
            e.hyperGammaMeanDN <- pars$gammaMeanDN[j2]
            e.hyperBetaDN <- pars$betaDN[j2]
            e.bf <- bayes.factor.denovo1(x =  dnData[, j2],
                                    N = pars$nfamily[j2],
                                    mu =  mutData[, j2],
                                    gamma.mean = e.hyperGammaMeanDN,
                                    beta = e.hyperBetaDN)
            bfDN[, j2] <- e.bf$BF
            llk1DN[, j2] <- e.bf$llk1
            llk0DN[, j2] <- e.bf$llk0
        }
        bfAll <- bfAll*apply(bfDN, 1, prod)
        llk1All <- llk1All*apply(llk1DN, 1, prod)
        llk0All <- llk0All*apply(llk0DN, 1, prod)
    }

    if (length(pars$gammaMeanCC) == 0) {
        message("No parameters for case-control data;  therefore, these categories are not calculated in this step.\n")
    } else {

        bfCC <- llk1CC <- llk0CC <- matrix(1, ncol = dim(caseData)[2], nrow = dim(caseData)[1])
        for (cc3 in 1:dim(bfCC)[2]){
            e.hyperGammaMeanCC <- pars$gammaMeanCC[cc3]
            e.hyperBetaCC <- pars$betaCC[cc3]
            e.nu <- 200
            t.case <- caseData[, cc3]
            t.control <- controlData[, cc3]
            e.rho <- e.nu*mean(t.case + t.control)/(pars$ncase[cc3] + pars$ncontrol[cc3])
            e.bf <- BayesFactorCC3(Nsample = list(ca = pars$ncase[cc3], cn = pars$ncontrol[cc3]),
                                   x.case = t.case, x.control = t.control,
                                   gamma.meanCC = e.hyperGammaMeanCC, betaCC = e.hyperBetaCC,
                                   rhoCC = e.rho, nuCC = e.nu)
            bfCC[, cc3] <- e.bf$BF
            llk1CC[, cc3] <- e.bf$llk1
            llk0CC[, cc3] <- e.bf$llk0
        }

        bfAll <- bfAll*apply(bfCC, 1, prod)
        llk1All <- llk1All*apply(llk1CC, 1, prod)
        llk0All <- llk0All*apply(llk0CC, 1, prod)

    }
  #  outData$pi0 <- pi0
 #   outData$llk1 <- llk1All
#    outData$
    outData$BF <- bfAll
    outData$PP <- bfAll*pi0/(1 - pi0 + pi0*bfAll)

    llk = sum(log(pi0*llk1All + (1 - pi0)*llk0All))
    

#    outData <- outData[order(-outData$BF),]
#    outData$qvalue <- Bayesian.FDR(outData$BF, 1 - pi0)$FDR

    return(list(dataFDR = outData, loglk = llk))
}

