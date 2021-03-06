DECO <- function(inputData,
                  Ndn = NULL,
                  Ncase = NULL, Ncontrol = NULL,
                  useBF = FALSE, # Using directly bayes factors
                  geneSetTest = FALSE, 
                  pi0 = NULL, ##If useBF=TRUE then pi0 has to be input.
                  nIteration = NULL,
                  ngsIteration = 20,
                  eLlk = 0.01, 
                  nThin = NULL, nCore = 1, nChain = 1,
                  nIteration2 = NULL,               
                  nThin2 = NULL,
                  estimationMethodforUsingBF = c("optimizing", "mcmc"),

                  hyperBetaDN0 = NULL,
                  hyperBetaCC0 = NULL,
                  hyper2GammaMeanDN = NULL, hyper2BetaDN = NULL, ##Priors for mean RRs: meanRR ~ gamma(hyper2GammaMeanDN, hyper2Be
                  hyper2GammaMeanCC = NULL, hyper2BetaCC = NULL,
                  alpha01SD = 2,
                  alpha01U = -1,
                  rho0 = NULL, nu0 = NULL,
                  upperPi0 = 0.5, lowerPi0 = 0, lowerBeta = 0, ##Lower and upper limits of pi: should be default
                  lowerHyperGamma = 1, lowerGamma = 1, #Should be default
                  betaPars = c(6.7771073, -1.7950864, -0.2168248), #Adjust beta's values: should be default
                  adjustHyperBeta = as.integer(1), ##1 if want to adjust beta, 0 if don't want to adjust beta
                  autoAdjustHyperBeta = FALSE,
                  drawHeatMap = FALSE, 
                  writeResult = FALSE,
                  resultDir = NULL,
                  initList = NULL,
                  seed0 = sample.int(.Machine$integer.max, 1))
{


###MCMC process
    geneName <- data.frame(inputData[, "Gene", drop = FALSE])
    dataDN <- data.frame(inputData[, grep("dn_", colnames(inputData)), drop = FALSE])
                                        #         colnames(dataDN) <- paste0("dn_", 1:dim(dataDN)[2])
    mutRate <- data.frame(inputData[, grep("mut_", colnames(inputData)), drop = FALSE])
                                        #        colnames(mutRate) <- paste0("mut_", 1:dim(mutRate)[2])
    dataCCcase <- data.frame(inputData[, grep("cc_case", colnames(inputData)), drop = FALSE])
                                        #       colnames(dataCCcase) <- paste0("cc_case", 1:dim(dataCCcase)[2])
    dataCCcontrol <- data.frame(inputData[, grep("cc_control", colnames(inputData)), drop = FALSE])

    dataBF <- data.frame(inputData[, grep("BF", colnames(inputData)), drop = FALSE])

    geneSetData0 <- data.frame(inputData[, grep("GS", colnames(inputData)), drop = FALSE])
    if (dim(geneSetData0)[2] != 0){
        geneSetData <- data.frame(Gene = geneName, GS = geneSetData0[, 1])
        geneSetTest <- TRUE
        
    }

    genesetInfo <- NULL ##Result of gene-set analysis
    mcmcData = NULL ##MCMC data of genetic parameters
    e.qOut <- NULL ##Output of q1 and q0
    pGS = NULL #p-value of the tested gene-set
    pH0.out = NULL ##probability of the null model
    pars0 = NULL ##Parameter results
    modelDN <- modelCC <- TRUE
    if (!is.null(resultDir))
        resultDir <- "."
    if (dim(dataCCcontrol)[2] == 0)
        dataCCcontrol = NULL
    if (dim(dataCCcase)[2] == 0){
        dataCCcase = NULL
        modelCC = FALSE #No case:control info
    }
    if (dim(dataDN)[2] == 0){
        dataDN = NULL
        modelDN = FALSE #No de novo info
    }
    if (dim(mutRate)[2] == 0)
        mutRate = NULL

    if (is.null(nIteration)){
        nIteration <- 5000
        message("No interation input; therefore, 5000 is used\n")
        
    }

    if (dim(dataBF)[2] != 0){
        useBF = TRUE
        message("Bayes Factors are used; therefore, it only takes some seconds to finish all analyses.\n")
        dataFDR <- inputData
        
    }
########Set a model: Both DN and CC, only DN or only CC
    if (modelDN & modelCC){
        modelName <- "DNandCCextTADA"
        message("\n=====Both de novo and case:control data are used=============\n")
    } else {
        if (modelDN){
            modelName <- "DNextTADA"
            message("\n=====Only de novo data are used=========\n")
            
        } else {
            modelName <- "CCextTADA"
            message("\n=====Only case:control data are used===============\n")                 
        }
        
    }
    modelName <- eval(parse(text = modelName))
###################
    
    if (is.null(nThin))
        nThin <- floor(nIteration/1000)
                                        #                  message("print(head(dataDN)) :", head(dataDN))

    ##There are two sections: using BFs (faster) or not using BFs (slow, using MCMC)
    if (!(useBF)){

        if (is.null(pi0)){
            alpha01Mean <- -2
        }
        message("\n=======MCMC sampling============\n")
        mcmcData <- extTADAmcmc(modelName = modelName,
                                        #                               geneSet = geneSet,
                                dataDN = dataDN,
                                mutRate = mutRate, Ndn = Ndn,
                                dataCCcase = dataCCcase, dataCCcontrol = dataCCcontrol, Ncase = Ncase, Ncontrol = Ncontrol,
                                nIteration = nIteration, nIteration2 = nIteration2,
                                alpha01Mean = alpha01Mean, alpha01SD = alpha01SD,
                                alpha01U = alpha01U,
                                nThin = nThin, nThin2 = nThin2, nCore = nCore, nChain = nChain,
                                hyperBetaDN0 = hyperBetaDN0, hyperBetaCC0 = hyperBetaCC0,
                                hyper2GammaMeanDN = hyper2GammaMeanDN, hyper2BetaDN = hyper2BetaDN, ##Priors for mean RRs: meanRR ~ gamma(hyper2GammaMeanDN, hyper2BetaDN)
                                hyper2GammaMeanCC = hyper2GammaMeanCC, hyper2BetaCC = hyper2BetaCC,
                                upperPi0 = upperPi0, lowerPi0 = lowerPi0, lowerBeta = lowerBeta, ##Lower and upper limits of pi: should be default
                                lowerHyperGamma = lowerHyperGamma, lowerGamma = lowerGamma, #Should be default
                                betaPars = betaPars, #Adjust beta's values: should be default
                                adjustHyperBeta = adjustHyperBeta, ##1 if want to adjust beta, 0 if don't want to adjust beta
                                autoAdjustHyperBeta =  autoAdjustHyperBeta)



###############Estimate genetic parameters
        message("\nEstimate genetic parameters from MCMC results")
############Add parameters
        a_parNames <- colnames(as.data.frame(mcmcData))
        parNames <- a_parNames[grep('alpha0', a_parNames)]

        if (!is.null(dataCCcontrol)){
            parNames <- c(parNames, paste0('hyperGammaMeanCC[', 1:dim(dataCCcontrol)[2], ']'))
            parNames <- c(parNames, paste0('hyperBetaCC[', 1:dim(dataCCcontrol)[2], ']'))
        }
        if (!is.null(dataDN)){
            parNames <- c(parNames, paste0('hyperGammaMeanDN[', 1:dim(dataDN)[2], ']'))
            parNames <- c(parNames, paste0('hyperBetaDN[', 1:dim(dataDN)[2], ']'))
        }
        

        message("Estimate parameters from MCMC results\n")
        
        pars0 <- estimatePars(pars = parNames,
                              mcmcResult = mcmcData)
        print(pars0)
        pars1 <- as.numeric(pars0[, 1])
        names(pars1) <- rownames(pars0)

        parsFDR <- list(alpha0 = as.numeric(pars1[grep("alpha0", names(pars1))]),
                        gammaMeanDN = as.numeric(pars1[grep("hyperGammaMeanDN", names(pars1))]),
                        betaDN = as.numeric(pars1[grep("hyperBetaDN", names(pars1))]),
                        gammaMeanCC = as.numeric(pars1[grep("hyperGammaMeanCC", names(pars1))]),
                        betaCC = as.numeric(pars1[grep("hyperBetaCC", names(pars1))]),

                        nfamily = Ndn,
                        ncase = Ncase,
                        ncontrol = Ncontrol
                        )

        message("\nCalculate posterior probabilities and FDRs")
        colnames(geneName) <- "Gene"
        dataOut <- calculateFDR(pars = parsFDR, #geneSet = geneSet,
                                dnData = dataDN, mutData = mutRate,
                                caseData = dataCCcase, controlData = dataCCcontrol,
                                geneName = geneName)

        dataFDR <- dataOut$dataFDR
        dataFDR[, 'FDR'] <- 1 - dataFDR[, 'PP']
        dataFDR <- dataFDR[order(dataFDR[, 'FDR']), ]
        fdr.x0 <- fdr.x1 <- dataFDR$FDR

        for (i in 1:length(fdr.x0))
            fdr.x1[i] <- sum(fdr.x0[1:i])/i
        dataFDR[, 'FDR'] <- fdr.x1

        loglk <- dataOut$loglk
        message("\nDraw heatmaps")
        outTime <-  format(Sys.time(), "%a_%b_%d_%H_%M_%S_%Y")
        if (drawHeatMap) {
            pdf(paste0(resultDir, "/heatMap", outTime, ".pdf"))
            allHyperGamma <- rownames(pars0[grep("hyperGammaMean", rownames(pars0)), ])
            for (i1 in 1:length(allHyperGamma))
                plotParHeatmap(pars = c("alpha0[1]", allHyperGamma[i1]), mcmcResult = mcmcData)
            dev.off()
        }
        if (writeResult){
            write.table(dataFDR, paste0(resultDir, "/Result_extTADA_PosteriorAndFDR", outTime, ".txt"),
                        row.names = FALSE, quote = FALSE)
            write.table(pars0,   paste0("Result_extTADA_estimatedPars", outTime, ".txt"), quote = FALSE)
        }
########################
        message(paste0("\nThe analysis is completed.\nIf you want to analyse steps seperately, please take a look at the example in the manual"))


        piTemp <- as.numeric(pars0[1, 1:3])
        piTemp <- exp(piTemp)/(1 + exp(piTemp))
        pi0 <- piTemp[1]
        piTemp <- c('pi0', piTemp, 0, 0, 0)

        print(piTemp)
        pars0 <- cbind(rownames(pars0), pars0)
        colnames(pars0)[1] <- c("Parameter")

        
        pars0 <- t(data.frame(pi0 = piTemp, t(pars0)))

        if (dim(geneSetData0)[2] != 0){
            print("Adding gene-set info")
            print(head(geneSetData, 2))
            print(head(dataFDR, 2))
            
            colnames(dataFDR)[1] <- c("Gene")
            dataFDR <- merge(dataFDR, geneSetData, by = 'Gene')
            print(head(dataFDR, 2))
            print("Done")
        }

        dataBF <- data.frame(dataFDR[, grep("BF", colnames(dataFDR)), drop = FALSE])
        useBF <- TRUE

    }


        
        pH0 <- p0.temp.dn <- p0.temp.cc <- rep(1, dim(inputData)[1])
        if (!is.null(dataDN)){
            p0.temp.dn <- pH0.model.DN(data.dn = dataDN, data.mut = mutRate, n.dn = Ndn)
            pH0 <- pH0*apply(p0.temp.dn, 1, prod)
        }
        if (!is.null(dataCCcase)){
            data.cc <- data.frame(dataCCcase, dataCCcontrol)
            p0.temp.cc <- pH0.model.CC(data.cc = data.cc,
                                       n.cc  = list(ncase = Ncase, ncontrol = Ncontrol),
                                       rho0 = rho0, nu0 = nu0)
            pH0 <- pH0*apply(p0.temp.cc, 1, prod)
        }
        pH0.out = data.frame(Gene = geneName, pH0 = pH0)
        Ngene = dim(inputData)[1]

###Optimizing function

        if (geneSetTest){
        message("\n========Bayes Factors are used in the estimation of gene-set information==========\n")

        if (is.null(pi0))
            stop("pi0 has to be input if BF is used")
        estimationMethodforUsingBF <- match.arg(estimationMethodforUsingBF)

        
    
        if (dim(dataFDR[, grep("FDR", colnames(dataFDR)), drop = FALSE])[2] != 0)
            dataFDR[, 'FDR0'] <- dataFDR[, 'FDR']
        if (dim(dataFDR[, grep("PP", colnames(dataFDR)), drop = FALSE])[2] != 0)
            dataFDR[, 'PP0'] <- dataFDR[, 'PP']

        bf0 = dataBF[, 1]
        gs0 = dataFDR[, 'GS'] #geneSet[, 1]

        logLK <- NULL
        e.qList <- NULL
        xGS <- dataFDR[, 'GS', drop = FALSE]
        piM <- c(1 - pi0, pi0)
        q00 <- apply(xGS, 2, function(x) length(x[x > 0])/length(x))
        e.q <- matrix(rep(q00, length(piM)), nrow = length(piM), byrow = TRUE) #Nmodel x Ngs
        e.qList <- logLK <- NULL

        for (kk1 in 1:ngsIteration){
#####################################################
                                        #E step
#####################################################
            E1 <- t(apply(cbind(bf0, xGS), 1, function(y){
                ##zij = 
                zVector <- c(piM[1], piM[2]*y[1])
                for (jz in 1:dim(xGS)[2]){
                    for (iz in 1:length(piM)){
                        zVector[iz] <- exp(sum(log(zVector[iz]) + log(dbinom(y[1 + jz], 1, prob = e.q[iz, jz]))))
                    }                }
                
                return(zVector)
            }))
            head(E1, 2)
            E2 <- t(apply(E1, 1, function(x) x/sum(x)))
            head(E2, 2)
###################################################
                                        # M step
###################################################
            for (iz in 2:length(piM)){
                for (jz in 1:dim(xGS)[2])
                    e.q[iz, jz] <- sum(E2[, iz]*xGS[, jz])/sum(E2[, iz])
            }
            e.qList[[kk1]] <- e.q
                                        # Log likelihood
            LK <- rowSums(E1)
            LK <- pH0*LK
            logLK[kk1] <- sum(log(LK))

            if (kk1 > 1)
                if (abs(logLK[kk1] - logLK[kk1 -1 ]) < eLlk)
                    break
        }
        
##################
#####################
        l0 <- logLK[1] 
        l1 <- logLK[length(logLK)]
        lChisq <- -2*(l0 - l1)
        pGS <- pchisq(abs(lChisq), df = 1, lower.tail = FALSE)
        genesetInfo = list(pGS = pGS, q01 = e.q, q00 = q00)
        
        dataFDR$BF <- bf0
        dataFDR$PP <- E2[, 2] #

        dataFDR$FDR <- 1 - dataFDR$PP
        dataFDR <- dataFDR[order(dataFDR$FDR), ]
        fdr.x0 <- dataFDR$FDR

        for (i in 1:length(fdr.x0))
            dataFDR$FDR[i] <- sum(fdr.x0[1:i])/i

##        q01 = 2
##        se.mle = 1
##        pars0 = rbind(c(pi0, pi0, pi0, 0, 0, 0), #c(alpha01Mean, alpha01Mean, alpha01Mean, 0, 0, 0),
##                      c(q01, q01 - 1.96*se.mle, q01 + 1.96*se.mle, se.mle, q01/se.mle, pGS))
##        rownames(pars0) <- c('pi0', 'p1GeneSet')
##        colnames(pars0) <- c("Mode", "lCI", "uCI", "pSD", "zScore",                             "pValue")

##        pars0 <- cbind(rownames(pars0), pars0)
##        colnames(pars0)[1] <- c("Parameter")

        loglk <- logLK
        e.qOut <- do.call(cbind, e.qList)

    }
    
    
    return(list(dataPP = dataFDR, pars = pars0, genesetInfo = genesetInfo, mcmcData = mcmcData, loglk = loglk, e.q = e.qOut, pH0 = pH0.out))


}

