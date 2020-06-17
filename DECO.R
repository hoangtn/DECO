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

#################################################################
# Some useful functions
#################################################################

# Fisher's method of combining p values
# Input: a vector of p-values
combpval <- function(p) {
  k <- length(p)
  T <- -2 * sum(log(p))
  return (1 - pchisq(T,2*k))
}

# Genomic control
# Input: a vector of p-values, the quantile (0.5 or 0.75)
genom.ctrl <- function(p, quant) {
  chisq.obs <- qchisq(1-p, 1) # convert p-values to chi-square statistics (dof 1)
  chisq.obs.quant <- quantile(chisq.obs, probs=quant, names=FALSE) # 0.75 if first quantile
  chisq.exp.quant <- qchisq(quant, 1)
  return (chisq.obs.quant / chisq.exp.quant)
}

# Bayesian FDR control (PMID:19822692, Section2.3)
# BF: a sorted vector of BFs (in decreasing order)
# pi0: the prior probability that the null model is true
# alpha: the FDR target
# Return: the q-value of each BF, and the number of findings with q below alpha. 
Bayesian.FDR <- function(BF, pi0, alpha=0.05) {
  # convert BFs to PPA (posterior probability of alternative model)
  pi <- 1-pi0
  q <- pi*BF/(1-pi+pi*BF) # PPA
  q0 <- 1 - q # posterior probability of null model
  
  # the FDR at each PPA cutoff
  n <- length(BF)
  FDR <- numeric(n)
  for (i in 1:n) FDR[i] <- sum(q0[1:i]) / i 
  
  # the cutoff
  t <- 1
  while (t <= length(q0) & mean(q0[1:t]) <= alpha) { t <- t+1 }
  return (list(FDR=FDR, ND=t))
}

# draw QQ plot
# Input: a vector of p-values
plotQQ <- function(p.obs) {
  obs <- -log10(p.obs)
  theo <- -log10(ppoints(length(obs)))
  qqplot(theo, obs, xlab=expression("Theoretical " * -log[10](p)), ylab=expression("Observed "*-log[10](p)))
  abline(0,1,col='red')
}

# Similar, but show QQ plot in the original scale (not log. scale)
plotQQ.unif <- function(p.obs) {
  obs <- (p.obs)
  theo <- (ppoints(length(obs)))
  qqplot(theo, obs, xlab="Theoretical p-values", ylab="Observed p-values")
  abline(0,1,col='red')
}

#################################################################
# Bayes Factor Computation for a Single Gene
#################################################################

# model evidence of de novo data: P(x_d|H_0) 
# Input: the count data x, the sample size N, the mutation rate mu 
evidence.null.dn <- function(x, N, mu) {
  return (dpois(x, 2*N*mu))
}

# model evidence of de novo data: P(x_d|H_1) 
# Input: the count data x, the sample size N, the mutation rate mu and the parameters
# Prior distribution of RR: gamma ~ Gamma(gamma.mean*beta, beta)
evidence.alt.dn <- function(x, N, mu, gamma.mean, beta) {
  return (dnbinom(x, gamma.mean*beta, beta/(beta+2*N*mu)))
}


# model evidence of case-control data: P(x_1,x_0|H_0) 
# Input: the count data x, the sample size N and the parameters
# Prior distribution of q|H0: Gamma(rho0, nu0)
evidence.null.CC <- function(x, N, rho0, nu0) {
  marglik0.ctrl.log <- log(dnbinom(x$cn, rho0, nu0/(nu0+N$cn)))
  marglik0.case.log <- log(dnbinom(x$ca, rho0+x$cn, (nu0+N$cn)/(nu0+N$cn+N$ca)))
  marglik0.log <- marglik0.ctrl.log + marglik0.case.log
  
  return (list(ctrl=exp(marglik0.ctrl.log), case=exp(marglik0.case.log), total=exp(marglik0.log)))
}

# model evidence of case-control data: P(x_1,x_0|H_1) 
# Input: the count data x, the sample size N and the parameters
# Prior distribution of RR: gamma ~ Gamma(gamma.mean*beta, beta)
# Prior distribution of q|H1: Gamma(rho1, nu1)
evidence.alt.CC <- function(x, N, gamma.mean, beta, rho1, nu1, q.lower=1e-8, q.upper=0.1, debug=FALSE) {
  integrand <- function(u) {
    q <- exp(u)
    return (dnbinom(x$ca, gamma.mean*beta, beta/(beta+N$ca*q)) * dgamma(q, rho1+x$cn, nu1+N$cn) * exp(u))
  }
  
  marglik1.ctrl <- dnbinom(x$cn, rho1, nu1/(nu1+N$cn))
  marglik1.case <- integrate(integrand, lower=log(q.lower), upper=log(q.upper))$value
  
  marglik1 <- marglik1.ctrl * marglik1.case
  
  #   return (exp(marglik1.ctrl.log+marglik1.case.log))
  return (list(ctrl=marglik1.ctrl, case=marglik1.case, total=marglik1))
}

# Bayes factor of the case-control data
# BF.cn and BF.ca: contribution from control and case data, respectively
# Prior distribution of RR: gamma ~ Gamma(gamma.mean*beta, beta)
# Prior distribution of q|H1: Gamma(rho1, nu1)
# Prior distribution of q|H0: Gamma(rho0, nu0)
bayes.factor.CC <- function(x, N, gamma.mean, beta, rho1, nu1, rho0, nu0) {
  marglik0.CC <- evidence.null.CC(x, N, rho0, nu0)
  marglik1.CC <- evidence.alt.CC(x, N, gamma.mean, beta, rho1, nu1)
  
  BF.cn <- marglik1.CC$ctrl / marglik0.CC$ctrl
  BF.ca <- marglik1.CC$case / marglik0.CC$case
  BF <- BF.cn * BF.ca
  
  return (list(BF=BF, BF.cn=BF.cn, BF.ca=BF.ca))
}

# Bayes factor of de novo counts of a gene 
# x: the de novo count
# N: the sample size (number of families)
# mu: the mutation rate (of this type of mutational events)
# Prior distribution of RR: gamma ~ Gamma(gamma.mean*beta, beta)
bayes.factor.denovo <- function(x, N, mu, gamma.mean, beta) {
  marg.lik0 <- dpois(x, 2*N*mu)
  marg.lik1 <- dnbinom(x, gamma.mean*beta, beta/(beta+2*N*mu))
  BF <- marg.lik1/marg.lik0
  
  return (BF)
}

# Bayes factor of the gene combining de novo and case-control
# x: a list of (dn, ca, cn), counts in de novo, cases and controls
# N: a list of (dn, ca, cn), sample sizes
# hyperpar: (gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0)
# Prior distribution of RR in de novo: gamma.dn ~ Gamma(gamma.mean.dn*beta.dn, beta.dn)
# Prior distribution of RR in C/C data: gamma.CC ~ Gamma(gamma.mean.CC*beta.CC, beta.CC)
# Prior distribution of q|H1: Gamma(rho1, nu1)
# Prior distribution of q|H0: Gamma(rho0, nu0)
bayes.factor <- function(x, N, mu, hyperpar, debug=FALSE) {
  gamma.mean.dn <- hyperpar[1]
  beta.dn <- hyperpar[2]
  gamma.mean.CC <- hyperpar[3]
  beta.CC <- hyperpar[4]
  rho1 <- hyperpar[5]
  nu1 <- hyperpar[6]
  rho0 <- hyperpar[7]
  nu0 <- hyperpar[8]
  
  BF.dn <- bayes.factor.denovo(x$dn, N$dn, mu, gamma.mean.dn, beta.dn)
  x.CC <- list(ca=x$ca, cn=x$cn)
  N.CC <- list(ca=N$ca, cn=N$cn)
  BF.CC <- bayes.factor.CC(x.CC, N.CC, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0)$BF
  BF <- BF.dn * BF.CC
  
  if (debug) {
    cat("BF.dn = ", BF.dn, "\n")
    cat("BF.CC = ", BF.CC, "\n")
  }
  
  return (BF)
}

#################################################################
# TADA-Denovo: analysis of de novo data 
#################################################################

# Genome-wide application of TADA for denovo data
# Input: counts, N, mu, gamma.mean, beta 
# counts: m x K matrix, where m is the number of gene, and K is the number of mutational categories. counts[i,j] is the number of de novo mutation in the j-th category of the i-th gene. 
# mu: mutation rates of the genes
# mu.frac: the fraction of each category (vector)
# gamma.mean, beta: the parameters of RR. Vector (one value for each category)
# Output: the statistics (BFs), an mxK matrix, and the total BFs (m-dim. vector)
TADA.denovo <- function(counts, N, mu, mu.frac, gamma.mean, beta) {
  m <- length(mu)
  K <- length(mu.frac)
  BF <- array(1, dim=c(m,K))
  BF.total <- numeric(m)
  
  # Compute BFs
  for (i in 1:m) {
    for (j in 1:K)  BF[i,j] <- bayes.factor.denovo(counts[i,j], N, mu[i]*mu.frac[j], gamma.mean[j], beta[j])
  }
  
  # Total BF per gene
  BF.total <- exp(rowSums(log(BF)))
  
  return (list(BF=BF, BF.total=BF.total))
}

# Compute permutation BFs of one gene: de novo only
# mu.gene: the mutation rates of a gene (K-dim. vector)
# N: the sample size
# l: the number of permutations
# gamma.mean, beta: RR of de novo mutations (vectors)
# Output: BF - l BFs from permutation; sample - permutate data
permute.gene.denovo <- function(mu.gene, N, l, gamma.mean, beta, dn.max=5) {
  K <- length(mu.gene)
  BF.gene <- array(1, dim=c(l,K))
  BF.gene.total <- numeric(l)
  
  # permutation of l times
  count.gene <- array(0, dim=c(l,K))
  for (j in 1:K) {
    # pre-compute the de novo table for the j-th category
    table.dn <- numeric(dn.max+1)
    for (k in 0:dn.max) {
      table.dn[k+1] <- bayes.factor.denovo(k, N, mu.gene[j], gamma.mean[j], beta[j])
    }
    
    # permutation
    count.gene[,j] <- rpois(l, 2*N*mu.gene[j])
    for (i in 1:l) {
      x <- count.gene[i,j]
      cond.range.dn <- (x <= dn.max)
      if (cond.range.dn==TRUE) {
        BF.gene[i,j] <- table.dn[x+1]
      } else {
        BF.gene[i,j] <- bayes.factor.denovo(x, N, mu.gene[j], gamma.mean[j], beta[j]) 
      }
    }
  }
  
  BF.gene.total <- exp(rowSums(log(BF.gene)))
  
  return (list(BF=BF.gene.total, sample=count.gene))
}

# Genome-wide application of TADA for denovo data: the difference with TADA.denovo is the report of p-values. 
# Input: counts, N, mu, mu.frac, gamma.mena, beta 
# counts: m x K matrix, where m is the number of gene, and K is the number of mutational categories. counts[i,j] is the number of de novo mutation in the j-th category of the i-th gene. 
# mu: mutation rates of the genes
# mu.frac: the fraction of each category (vector)
# gamma.mean, beta: the parameters of RR. Vector (one value for each category)
# l: the number of permutations to obtain the null distribution
# Output: the statistics (BFs), an mxK matrix, and the total BFs (m-dim. vector). pval: the p-values of all genes. BF.null: the null distribution of BFs. 
TADAp.denovo <- function(counts, N, mu, mu.frac, gamma.mean, beta, l=100, dn.max=5) {
  m <- length(mu)
  K <- length(mu.frac)
  BF <- array(1, dim=c(m,K))
  BF.total <- numeric(m)
  pval <- numeric(m)
  
  # Compute BFs
  rs <- TADA.denovo(counts, N, mu, mu.frac, gamma.mean, beta)
  BF <- rs$BF
  BF.total <- rs$BF.total
  
  # Create the null distribution
  BF.null <- numeric(m*l)
  for (i in 1:m) {
    BF.null[((i-1)*l+1):(i*l)] <- permute.gene.denovo(mu[i]*mu.frac, N, l, gamma.mean, beta, dn.max=dn.max)$BF
  }
  
  # p-values of each gene
  BF.null.srt <- sort(BF.null, decreasing=TRUE)
  pval <- findInterval(-BF.total, -BF.null.srt)/length(BF.null.srt)
  pval[pval==0] <- 0.5/length(BF.null.srt)
  
  return (list(BF=BF, BF.total=BF.total, pval=pval, BF.null=BF.null))
}

#################################################################
# TADA: analysis of de novo and inherited data
#################################################################

# Genome-wide application of TADA 
# counts: m x 3K matrix, where m is the number of gene, and K is the number of mutational categories. Each category has three numbers: de novo, case and control. 
# N: sample size, three values for de novo, case and control
# mu: mutation rates of the genes
# mu.frac: the fraction of each category (vector)
# hyperpar: 8*K matrix, where each row is (gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0), and each column corresponds to one mutation type
# denovo.only: whether using only de novo data (Boolean vector)
# pi.gene: for each gene, the estimated fractions of causal variants, one for each class of variants. These fractions will be used to set gene-specific RR (case-control)
# Output: the statistics (BFs), an mxK matrix, and the total BFs (m-dim. vector)
TADA <- function(counts, N, mu, mu.frac, hyperpar, denovo.only=FALSE, pi.gene=1) {
  m <- length(mu)
  K <- length(mu.frac)
  BF <- array(1, dim=c(m,K))
  BF.total <- numeric(m)
  
  if (length(denovo.only)==1) { denovo.only <- rep(denovo.only, m) }
  if (length(pi.gene)==1) { pi.gene <- array(1, dim=c(m,K))}
  
  gamma.mean.dn <- hyperpar[1,]
  beta.dn <- hyperpar[2,]
  gamma.mean.CC <- hyperpar[3,]
    
  # Compute BFs: BF[i,j] is the BF of the i-th gene in the j-th category
  for (i in 1:m) {
    if (denovo.only[i]==FALSE) {
      for (j in 1:K)  {
        # set hyperparameters
        hyperpar.gene <- hyperpar[,j]
        RR.product <- hyperpar.gene[3]*hyperpar.gene[4]
        hyperpar.gene[3] <- hyperpar.gene[3]*pi.gene[i,j] + (1-pi.gene[i,j])
        hyperpar.gene[4] <- RR.product/hyperpar.gene[3]
        
        # compute BF  
        start <- 3*(j-1)+1
        x <- list(dn=counts[i, start], ca=counts[i, start+1], cn=counts[i, start+2])
        BF[i,j] <- bayes.factor(x, N, mu[i]*mu.frac[j], hyperpar.gene)
      }
    } else {
      for (j in 1:K) {
        start <- 3*(j-1)+1
        x <- counts[i,start]
        BF[i,j] <- bayes.factor.denovo(x, N$dn, mu[i]*mu.frac[j], gamma.mean.dn[j], beta.dn[j])
      } 
    }
  }
  
  # Total BF per gene
  BF.total <- exp(rowSums(log(BF)))
  
  return (list(BF=BF, BF.total=BF.total))
}

# Compute permutation BFs of one gene
# mu.gene: the mutation rates of a gene (K-dim. vector), and the case-control counts (to be permuted): vectors (one value per category)
# N: sample size, three values for de novo, case and control
# l: number of permutations
# gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0: parameters. 
# table.CC: precomputed BFs. A 3-dim table, table.CC[i, j, k] stores the BF of (j-1, k-1) in the i-th category. 
# Output: BF - l BFs from permutation; sample.dn, sample.ca, sample.cn - permutate data
permute.gene <- function(mu.gene, count.ca, count.cn, N, l, gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0, table.CC, dn.max=5) {
  K <- length(mu.gene)
  BF.dn <- array(1, dim=c(l,K))
  BF.CC <- array(1, dim=c(l,K))
  BF.gene.total <- numeric(l)
  ca.max <- dim(table.CC)[1] - 1
  cn.max <- dim(table.CC)[2] - 1
  
  # permutation l times for each category
  sample.dn <- array(0, dim=c(l,K))
  sample.ca <- array(0, dim=c(l,K))
  sample.cn <- array(0, dim=c(l,K))
  for (j in 1:K) {
    # pre-compute the de novo table for the j-th category
    table.dn <- numeric(dn.max+1)
    for (k in 0:dn.max) {
      table.dn[k+1] <- bayes.factor.denovo(k, N$dn, mu.gene[j], gamma.mean.dn[j], beta.dn[j])
    }
    
    # generate permutation data
    sample.dn[,j] <- rpois(l, 2*N$dn*mu.gene[j])
    sample.ca[,j] <- rhyper(l, count.ca[j] + count.cn[j], N$ca+N$cn-count.ca[j]-count.cn[j], N$ca)
    sample.cn[,j] <- count.ca[j] + count.cn[j] - sample.ca[,j]
    
    # compute the BFs
    for (i in 1:l) {
      x <- list(dn=sample.dn[i,j], ca=sample.ca[i,j], cn=sample.cn[i,j])
      cond.range.dn <- (x$dn <= dn.max)
      if (cond.range.dn==TRUE) {
        BF.dn[i,j] <- table.dn[x$dn+1]
      } else {
        BF.dn[i,j] <- bayes.factor.denovo(x$dn, N$dn, mu.gene[j], gamma.mean.dn[j], beta.dn[j]) 
      }
      
      cond.range.CC <- (x$ca <= ca.max) & (x$cn <= cn.max)
      if ( cond.range.CC==TRUE ) {
        BF.CC[i,j] <- table.CC[j, x$ca + 1, x$cn + 1]
      } else {
        BF.CC[i,j] <- bayes.factor.CC(x, N, gamma.mean.CC[j], beta.CC[j], rho1[j], nu1[j], rho0[j], nu0[j])$BF
      }
    }
  }
  
  BF.gene <- BF.dn * BF.CC
  BF.gene.total <- exp(rowSums(log(BF.gene)))
  
  return (list(BF=BF.gene.total, sample.dn=sample.dn, sample.ca=sample.ca, sample.cn=sample.cn))
}

# Genome-wide application of TADA 
# counts: m x 3K matrix, where m is the number of gene, and K is the number of mutational categories. Each category has three numbers: de novo, case and control. 
# N: sample size, three values for de novo, case and control
# mu: mutation rates of the genes
# mu.frac: the fraction of each category (vector)
# hyperpar: 8*K matrix, where each row is (gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0), and each column corresponds to one mutation type
# l: the number of permutations to obtain the null distribution
# dn.max, ca.max, cn.max: if counts are below these values, the BFs will be pre-computed. 
# Output: the statistics (BFs), an mxK matrix, and the total BFs (m-dim. vector). pval: the p-values of all genes. BF.null: the null distribution of BFs.
TADAp <- function(counts, N, mu, mu.frac, hyperpar, l=100, dn.max=5, ca.max=10, cn.max=10) {
  gamma.mean.dn <- hyperpar[1,]
  beta.dn <- hyperpar[2,]
  gamma.mean.CC <- hyperpar[3,]
  beta.CC <- hyperpar[4,]
  rho1 <- hyperpar[5,]
  nu1 <- hyperpar[6,]
  rho0 <- hyperpar[7,]
  nu0 <- hyperpar[8,]
  
  m <- length(mu)
  K <- length(mu.frac)
  BF <- array(1, dim=c(m,K))
  BF.total <- numeric(m)
  pval <- numeric(m)
  
  # Compute BFs of the observed data
  rs <- TADA(counts, N, mu, mu.frac, hyperpar)
  BF <- rs$BF
  BF.total <- rs$BF.total
  
  # Pre-compute the bayes-factors of the case-control data
  table.CC <- array(1, dim=c(K, (ca.max+1), (cn.max+1)))
  for (j in 1:K) {
    for (x1 in 0:ca.max) {
      for (x0 in 0:cn.max) {
        x <- list(ca=x1,cn=x0)
        table.CC[j, x1+1, x0+1] <- bayes.factor.CC(x, N, gamma.mean.CC[j], beta.CC[j], rho1[j], nu1[j], rho0[j], nu0[j])$BF
      }
    }
  }
  
  # Create the null distribution
  BF.null <- numeric(m*l)
  for (i in 1:m) {
#     print(i)
    BF.null[((i-1)*l+1):(i*l)] <- permute.gene(mu[i]*mu.frac, counts[i, 3*(1:K)-1], counts[i, 3*(1:K)], N, l, gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0, table.CC, dn.max=dn.max)$BF
  }
  
  # p-values of each gene
  BF.null.srt <- sort(BF.null, decreasing=TRUE)
  pval <- findInterval(-BF.total, -BF.null.srt)/length(BF.null.srt)
  pval[pval==0] <- 0.5/length(BF.null.srt)
  
  return (list(BF=BF, BF.total=BF.total, pval=pval, BF.null=BF.null))
}

#################################################################
# MOM estimation of hyperprior parameters from de novo data
#################################################################

# Prob. of having d or more de novo mutations under H1 
# Use simulation, but could also use analytic form 
multihit.prob <- function(N, mu, gamma.mean, beta, d=2, S=100) {
  p <- numeric(S)
  gamma <- rgamma(S, gamma.mean*beta, rate=beta)
  for (i in 1:S) {
    p[i] <- 1 - ppois(d-1, 2*N*mu*gamma[i])
  }
  return (mean(p))
}

# Estimate the number of multihit genes in a genome. 
# d: the parameter of the multiplicity test. 
# Returns: M0 - the number of multihit genes from non-risk genes; M1 - the number from risk genes. 
count.multihit <- function(N, mu, pi, gamma.mean, beta, d=c(2,3), S=2) {
  m <- length(mu)
  M0 <- numeric(length(d))
  M1 <- numeric(length(d))
  
  # M1: the number of causal genes having d or more de novo mutations
  p.alt <- array(0, dim=c(m, length(d)))  # p1[i,j] = P(X_i >= d_j|H1)
  for (i in 1:m) {
    for (j in 1:length(d)) {
      p.alt[i,j] <- multihit.prob(N, mu[i], gamma.mean, beta, d=d[j], S=S)
    }
  }
  for (j in 1:length(d)) { 
    M1[j] <- m * pi  * mean(p.alt[,j]) 
  }
  
  # M0: the number of non-causal genes having d or more de novo mutations
  p.null <- array(0, dim=c(m, length(d)))  # p0[i,j] = P(X_i >= d_j|H0)
  for (i in 1:m) {
    for (j in 1:length(d)) {
      p.null[i,j] <- 1 - ppois(d[j] - 1, 2*N*mu[i])
    }
  }
  for (j in 1:length(d)) { 
    M0[j] <- m * (1-pi) * mean(p.null[,j]) 
  }
  
  result <- data.frame(d=d, M0=M0, M1=M1)
  return (result)
}

# Estimating relative risk and the number of multiple hits from de novo data
# Input: sample size (N), mutation rates of all genes (mu), observed number of de novo events (C), beta (parameter of the prior distribution of gamma), k (number of disease genes)
# Output: the average relative risk (gamma.mean), the expected number of multi-hit genes (M)
denovo.MOM <- function(N, mu, C, beta, k) {
  m <- length(mu) # number of genes
  
  # enrichment of de novo events
  nu <- C / (2 * N * sum(mu))
  
  # MOM estimator of gamma.mean
  gamma.mean <- (nu-1)*m/k +1
  
  # expected M (choose d = 2)
  rs <- count.multihit(N, mu, k/m, gamma.mean, beta, d=2)    
  M <- sum(rs$M1) + sum(rs$M0)
  
  return (list(gamma.mean=gamma.mean, M=M))
}

#################################################################
# Empirical Bayes estimation of hyperprior parameters
#################################################################

# Evalaute the marginal log-likelihood at given parameters
# pi: the fraction of causal genes
# counts: the count data (of one mutational category), a date frame
# prior.weight: putting a prior so that q1 and q0 tend to be close (recommended value: 5000). Default 0. 
# Output: the negative log-likelihood, and posterior, BF of each gene
marginal <- function(hyperpar, pi, counts, N, mu, prior.weight=0, denovo.only=FALSE, debug=FALSE) {
  gamma.mean.dn <- hyperpar[1]
  beta.dn <- hyperpar[2]
  gamma.mean.CC <- hyperpar[3]
  beta.CC <- hyperpar[4]
  rho1 <- hyperpar[5]
  nu1 <- hyperpar[6]
  rho0 <- hyperpar[7]
  nu0 <- hyperpar[8]
  
  n <- nrow(counts)
  prob <- numeric(n)
  posterior <- numeric(n)
  bf <- numeric(n)
  for (i in 1:n)  {
    if (debug) cat("i = ", i, "\tdn = ", counts[i,]$dn, "\tca = ", counts[i,]$ca, "\tcn = ", counts[i,]$cn, "\n")
    if (denovo.only==TRUE) {
      prob.M1 <- evidence.alt.dn(counts[i,"dn"], N$dn, mu[i], gamma.mean.dn, beta.dn) 
      prob.M0 <- evidence.null.dn(counts[i,"dn"], N$dn, mu[i])
    } else {
      prob.M1 <-  evidence.alt(counts[i,],N,mu[i],gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1)
      prob.M0 <-  evidence.null(counts[i,],N,mu[i],rho0, nu0) 
    }
    if (debug) cat("prob.M1 = ", prob.M1,"\n")
    if (debug) cat("Prob.M0 = ", prob.M0, "\n")
    prob[i] <- pi * prob.M1 + (1 - pi) * prob.M0
    posterior[i] <- pi*prob.M1 / (prob[i])
    bf[i] <- prob.M1 / prob.M0
  }
  marg.nll <- -sum(log(prob[!is.na(prob)]))
  
  # add the prior
  if (prior.weight > 0) {
    u <- prior.weight
    q1.mean <- rho1/nu1
    q0.mean <- rho0/nu0
    marg.nll <- marg.nll - log(dgamma(q1.mean, q0.mean*u, u)) 
  }
  
  result <- list(prob=prob, marginal=marg.nll, posterior=posterior, bf=bf)
  return (result)
}

# Convert the full set of parameters (pi, hyperpar) to a subset based on the option
# est.option: a Boolean vector (size 8), one for each hyperpar. If FALSE, the corresponding parameter will be fixed. 
# est.pi: Boolean. 
fullpar2subpar <- function(hyperpar, pi, est.option, est.pi) {
  subpar <- NULL
  for (i in 1:length(est.option)) {
    if (est.option[i]) subpar <- c(subpar, hyperpar[i])
  }
  
  if (est.pi==TRUE)  subpar <- c(subpar, pi)
  
  return (subpar)
}

# Convert the subset of parameters (subpar) to the full set based on the option and intial values
# hyperpar.init, pi.init: used to set parameters not to be estimated
subpar2fullpar <- function(subpar, hyperpar.init, pi.init, est.option, est.pi) {
  hyperpar <- numeric(length(hyperpar.init))
  curr <- 1
  for (i in 1:length(est.option)) {
    if (est.option[i]) { hyperpar[i] <- subpar[curr]; curr <- curr+1 }
    else hyperpar[i] <- hyperpar.init[i]
  }
  if (est.pi==TRUE)  pi <- subpar[curr]
  else pi <- pi.init
  
  return (list(hyperpar=hyperpar, pi=pi))
}

# Empirical Bayes estimation of hyperparameters 
# counts: the count data - the number of mutations in cases, controls and de novo data respectively. 
# N: sample size (dn, ca, cn, respectively). 
# mu: vector of mutation rates
# lower, upper: the search range of the parameters 
# est.option: a Boolean vector specifying whether to estimate each parameter. If FALSE, use the corresponding parameter in hyperpar.init
# est.pi: whether to estimate pi. If FALSE, use the given value of pi; if TRUE, use the given value as initial
# prior.weight: putting a prior so that q1 and q0 tend to be close (recommended value: 5000). Default 0. 
# Output: the hyperparameters and the NLL (negative log likelihood) at the parameters. 
empBayes <- function(counts, N, mu, hyperpar.init, pi.init, lower, upper, lower.pi=1e-10, upper.pi=1, est.option=rep(FALSE, 8), est.pi=FALSE, prior.weight=0, debug=FALSE) {  
  # marginal likelihood function (parameters in log-scale)
  marginal.loglike <- function(subpar.log) {
    subpar <- exp(subpar.log)
    allpar <- subpar2fullpar(subpar, hyperpar.init, pi.init, est.option, est.pi)
    hyperpar <- allpar$hyperpar
    pi <- allpar$pi
    result <- marginal(hyperpar,  pi, counts, N, mu, prior.weight=prior.weight)$marginal  
    
    if (debug) {
      cat("Parameters = (", paste(subpar, collapse=", "), ")\tValue = ", result, "\n", sep="") 
    }
    
    return (result)
  }
  
  # no parameter to estimate
  if ( sum(est.option==TRUE) == 0 & est.pi == FALSE ) {
    value <- marginal(hyperpar.init, pi.init, counts, N, mu, prior.weight=prior.weight)$marginal
    return (list(hyperpar=hyperpar.init, value=value))
  }
  
  # initialize the sub. parameters
  subpar.init <- fullpar2subpar(hyperpar.init, pi.init, est.option, est.pi)
  sublower <- fullpar2subpar(lower, lower.pi, est.option, est.pi)
  subupper <- fullpar2subpar(upper, upper.pi, est.option, est.pi)
  
  # maximization of marginal likelihood
  like.optim <- optim(log(subpar.init), marginal.loglike, method="L-BFGS-B", lower=log(sublower), upper=log(subupper))
  subpar.est.log <- like.optim$par
  value.max <- like.optim$value
  allpar.est <- subpar2fullpar(exp(subpar.est.log), hyperpar.init, pi.init, est.option, est.pi)
  hyperpar.est <- allpar.est$hyperpar
  pi.est <- allpar.est$pi
  
  if (est.pi==TRUE) { return (list(pi=pi.est, hyperpar=hyperpar.est, value=value.max)) }
  else { return (list(hyperpar=hyperpar.est, value=value.max)) }
}

# Empirical Bayes estimation of hyperparameters (gamma.mean, beta and pi) for de novo data
# counts: count data - the number of de novo mutations. 
# N - sample size (trios). 
# mu: vector of mutation rate
# lower, upper: the search range of the parameters 
empBayes.denovo <- function(counts, N, mu, gamma.mean.init, beta.init, pi.init, lower, upper, lower.pi=1e-10, upper.pi=1, debug=FALSE) {  
  # marginal likelihood function (parameters in log-scale)
  marginal.loglike <- function(par.log) {
    par <- exp(par.log)
    gamma.mean <- par[1]
    beta <- par[2]
    pi <- par[3]
    hyperpar <- c(gamma.mean, beta, 0, 0, 0, 0, 0, 0)
    m <- length(counts)
    counts.full <- data.frame(dn=counts, ca=numeric(m), cn=numeric(m))
    N.full <- list(dn=N, ca=0, cn=0)
    result <- marginal(hyperpar, pi, counts.full, N.full, mu, denovo.only=TRUE)$marginal  
    
    if (debug) {
      cat("Parameters = (", paste(par, collapse=", "), ")\tValue = ", result, "\n", sep="") 
    }
    
    return (result)
  }
  
  # initialize the parameters
  par.init <- c(gamma.mean.init, beta.init, pi.init)
  par.lower <- c(lower, lower.pi)
  par.upper <- c(upper, upper.pi)
  
  # maximization of marginal likelihood
  like.optim <- optim(log(par.init), marginal.loglike, method="L-BFGS-B", lower=log(par.lower), upper=log(par.upper))
  par.est <- exp(like.optim$par)
  value.max <- like.optim$value
  gamma.mean.est <- par.est[1]
  beta.est <- par.est[2]
  pi.est <- par.est[3]
  
  return (list(gamma.mean=gamma.mean.est, beta=beta.est, pi=pi.est, value=value.max)) 
}

#################################################################
# Functions for simulation
#################################################################

# Generate simulation data of a set of genes (multiple mutational categories): de novo mutations only
# N: sample size (number of trios)
# mu: mutation rate of each gene (a vector)
# mu.frac: for each type of mutation, its fraction of the gene's mutation rate
# pi: the fraction of risk genes 
# gamma.mean, beta: Relative risk of de novo mutation: gamma|M1 ~ Gamma(gamma.mean*beta, beta). Vectors (one per category) 
# Output: sample matrix (m by K), where m is the number of genes and K the number of variant categories. sample.info: more information of the samples, including the indicator (risk gene or not) and the RR. 
simulator.denovo <- function(N, mu, mu.frac, pi, gamma.mean, beta) {
  m <- length(mu) # number of genes
  K <- length(mu.frac) # number of mutational categories
  
  z <- rbinom(m, 1, pi)
  gamma <- array(1, dim=c(m,K))
  x <- array(0, dim=c(m,K))
  k <- sum(z==1)
  for (j in 1:K) {
    gamma[z==1, j] <- rgamma(k, gamma.mean[j]*beta[j], beta[j])
    x[,j] <- rpois(m, 2 * mu * mu.frac[j] * gamma[,j] * N)
  }
  
  sample.info <- cbind(mu, z, gamma, x)
  
  return (list(sample=x, sample.info=sample.info))
}

# Generate simulation data of a set of genes (multiple mutational categories)
# N: sample size (number of trios)
# mu: mutation rate of each gene (a vector)
# mu.frac: for each type of mutation, its fraction of the gene's mutation rate
# pi: the fraction of risk genes 
# gamma.mean.dn, beta.dn: Relative risk of de novo mutation: gamma|M1 ~ Gamma(gamma.mean.dn*beta.dn, beta.dn). Vectors.
# gamma.mean.CC, beta.CC: Relative risk of inherited mutation (case/control): gamma.CC|M1 ~ Gamma(gamma.mean.CC*beta.CC, beta.CC). Vectors
# Frequency parameter of risk genes: q|M1 ~ Gamma(rho1, nu1)
# Frequency parameter of non-risk genes: q|M0 ~ Gamma(rho0, nu0)
# tradeoff option: if TRUE, implement q-gamma tradeoff (i.e. higher gamma means lower q). Suppose, gamma_i is the RR, then q_i is proportional to mu_i / gamma_i, where the constant is determined from the mean of q, mu and gamma. 
# Output: sample matrix (m by 3K), where m is the number of genes and K the number of variant categories.  sample.info: more information of the samples, including the indicator (risk gene or not) and the RR.
simulator <- function(N, mu, mu.frac, pi, gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0, tradeoff=FALSE) {
  m <- length(mu) # number of genes
  K <- length(mu.frac) # number of mutational categories
  
  # the tradeoff parameter (delta:=mu.mean/q.mean)
  delta <- mean(mu) * mu.frac / (rho0 / nu0)
  
  z <- rbinom(m, 1, pi)
  gamma.dn <- array(1, dim=c(m,K))
  gamma.CC <- array(1, dim=c(m,K))
  q <- array(0, dim=c(m,K))
  x <- array(0, dim=c(m,3*K))
  k <- sum(z==1)
  for (j in 1:K) {
    # sample de novo 
    gamma.dn[z==1, j] <- rgamma(k, gamma.mean.dn[j]*beta.dn[j], beta.dn[j])
    col <- 3*(j-1)+1
    x[,col] <- rpois(m, 2 * mu * mu.frac[j] * gamma.dn[,j] * N$dn)
    
    # sample case-control
    gamma.CC[z==1, j] <- rgamma(k, gamma.mean.CC[j]*beta.CC[j], beta.CC[j])
    q[z==0, j] <- rgamma(m-k, rho0[j], nu0[j])
    if (tradeoff==FALSE) {
      q[z==1, j] <- rgamma(k, rho1[j], nu1[j])
    } else {
      q[z==1, j] <- mu[z==1] * mu.frac[j] / (delta[j] * gamma.CC[z==1, j])
    }
    x[,col+1] <- rpois(m, q[,j] * gamma.CC[,j] *N$ca)
    x[,col+2] <- rpois(m, q[,j] * N$cn)
    
  }
  
  sample.info <- cbind(mu, z, gamma.dn, gamma.CC, q, x)
  
  return (list(sample=x, sample.info=sample.info))
}

# Evaluation of the power (FDR) of TADA.denovo
# N: sample size, the number of families of the de novo study
# mu: mutation rate of each gene (a vector)
# mu.frac: for each type of mutation, its fraction of the gene's mutation rate
# pi: the fraction of risk genes 
# gamma.mean, beta: RR parameters for de novo mutations (vectors)
# gamma.mean.est, beta.est: the parameters used by TADA.denovo
# FDR: the FDR level to be controled
# Output: the number of discoveries at the specified FDR level
eval.TADA.denovo <- function(N, mu, mu.frac, pi, gamma.mean, beta, gamma.mean.est, best.est, FDR=0.1) {
  # sample the simulation data
  sample <- simulator.denovo(N, mu, mu.frac, pi, gamma.mean, beta)$sample
  
  # run TADA.denovo
  sampleBF <- TADA.denovo(sample, N, mu, mu.frac, gamma.mean.est, beta.est)$BF.total
  
  # Bayesian FDR control
  M1 <- Bayesian.FDR(sort(sampleBF, decreasing=TRUE), 1-pi, FDR)$ND
  
  return (M1)
}

# Evaluation of the power (FDR) of TADA
# N: sample size, the number of families of the de novo study
# mu: mutation rate of each gene (a vector)
# mu.frac: for each type of mutation, its fraction of the gene's mutation rate
# pi: the fraction of risk genes 
# gamma.mean.dn, beta.dn: RR parameters for de novo mutations (vectors)
# gamma.mean.CC, beta.CC: RR parameters for inherited mutations (vectors)
# Frequency parameter of risk genes: q|M1 ~ Gamma(rho1, nu1)
# Frequency parameter of non-risk genes: q|M0 ~ Gamma(rho0, nu0)
# hyperpar.est: the parameters used by TADA
# FDR: the FDR level to be controled
# Output: the number of discoveries at the specified FDR level
eval.TADA <- function(N, mu, mu.frac, pi, gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0, hyperpar.est, FDR=0.1, tradeoff=FALSE) {
  # sample the simulation data
  sample <- simulator(N, mu, mu.frac, pi, gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0, tradeoff=tradeoff)$sample
  
  # run TADA
  sampleBF <- TADA(sample, N, mu, mu.frac, hyperpar.est)$BF.total
  
  # Bayesian FDR control
  M1 <- Bayesian.FDR(sort(sampleBF, decreasing=TRUE), 1-pi, FDR)$ND
  
  # sample information (counts and BFs)
  sim.results <- data.frame(counts=sample, BF=sampleBF)
  sim.results <- sim.results[order(-sim.results$BF),]
  
  return (list(M1=M1, sample=sim.results))
}


bayes.factor.denovo1 <- function(x, N, mu, gamma.mean, beta) {
  marg.lik0 <- dpois(x, 2*N*mu)
  marg.lik1 <- dnbinom(x, gamma.mean*beta, beta/(beta+2*N*mu))
  BF <- marg.lik1/marg.lik0
  
  return (list(BF = BF, llk1 = marg.lik1, llk0 = marg.lik0))
}

BayesFactorCC3 <- function(x.case, x.control, Nsample,
                           gamma.meanCC, betaCC, rhoCC, nuCC){
    gAll <- range(rgamma(10000, gamma.meanCC*betaCC, rate = betaCC))
    gLower = gAll[1]; gUpper = gAll[2]

                                        #    print(range(gAll))
    altCC <- apply(cbind(x.case, x.control), 1, function(y){
        x2 <- list(ca = y[1], cn = y[2])
        evidence.alt.cc3 <- function(x = x2, N = Nsample, gamma.mean = gamma.meanCC, beta = betaCC,
                                     rho1 = rhoCC, nu1 = nuCC) {
            bControl <- log(dnbinom(x$cn, size = rho1, prob = nu1/(N$cn + nu1)))
                                        #                print(bControl)
            fCase <- function(gGamma) {
                dnbinom(x$ca, size = rho1 + x$cn, prob = (N$cn + nu1)/(N$cn + nu1 + N$ca*gGamma))*dgamma(gGamma, gamma.mean*betaCC, rate = betaCC)
            }
            bCase <- log(integrate(fCase, lower = gLower, upper = gUpper, stop.on.error = FALSE)$value)
                                        #print(bCase)
            return(exp(bCase + bControl))
        }
        t1 <- evidence.alt.cc3()

        return(t1)
    })
                                        #    print(altCC)
    nullCC <- apply(cbind(x.case, x.control), 1, function(y, rho1 = rhoCC, nu1 = nuCC, N = Nsample){
        x <- list(ca = y[1], cn = y[2])
        bControl <- log(dnbinom(x$cn, size = rho1, prob = nu1/(N$cn + nu1)))
        bCase <- log(dnbinom(x$ca, rho1 + x$cn, (N$cn + nu1)/(N$cn + nu1 + N$ca)))


        t1 <- exp(bCase + bControl)

        return(t1)
    })
                                        #    print(nullCC)

    tempBF <- altCC/nullCC #ifelse((x.case == 0) & (x.control == 0), 1, altCC/nullCC)
    return(list(BF = tempBF, llk1 = altCC, llk0 = nullCC))
}


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

##Calculate posterior probabilities
calculatePP <- function(pars,
                         caseData = NULL,
                         controlData = NULL,
                         dnData = NULL,
                         mutData = NULL,
                         geneName, geneSet){
    outData <- data.frame(geneName, geneSet)
    colnames(outData) <- c("GENE", paste0("GeneSet", 1:dim(geneSet)[2]))
    if (!is.null(dnData))
        outData <- cbind(outData, dnData)
    if (!is.null(mutData))
        outData <- cbind(outData, mutData)
    if (!is.null(caseData))
        outData <- cbind(outData, caseData)
    if (!is.null(controlData))
        outData <- cbind(outData, controlData)


    ##Calculate pi
    sumGeneSet = rep(pars$alphaValue[1], length(geneName))
    for (j in 1:dim(geneSet)[2])
        sumGeneSet = sumGeneSet + pars$alphaValue[j+1]*geneSet[, j]

    pi0 = exp(sumGeneSet)/(1 + exp(sumGeneSet))


    bfAll <- rep(1, dim(outData)[1])

    if ( length(pars$gammaMeanDN) == 0) {
        message("No parameters for de novo data; therefore, these categories are not calculated in this step.\n")
        }  else {
        bfDN <- matrix(1, nrow = dim(dnData)[1], ncol = dim(dnData)[2]) ##De novo bayes factors
        for (j2 in 1:dim(bfDN)[2]) {
            e.hyperGammaMeanDN <- pars$gammaMeanDN[j2]
            e.hyperBetaDN <- pars$betaDN[j2]
            e.bf <- bayes.factor.denovo(x =  dnData[, j2],
                                    N = pars$nfamily[j2],
                                    mu =  mutData[, j2],
                                    gamma.mean = e.hyperGammaMeanDN,
                                    beta = e.hyperBetaDN)
            bfDN[, j2] <- e.bf
        }
        bfAll <- bfAll*apply(bfDN, 1, prod)
    }

    if (length(pars$gammaMeanCC) == 0) {
        message("No parameters for case-control data;  therefore, these categories are not calculated in this step.\n")
    } else {

        bfCC <- matrix(1, ncol = dim(caseData)[2], nrow = dim(caseData)[1])
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
#            e.bf0 <- bayes.factor.CC(x

            
            bfCC[, cc3] <- e.bf
        }

        bfAll <- bfAll*apply(bfCC, 1, prod)
    }

    outData[, 'pi0'] <- pi0
    outData$BF <- bfAll
    outData[, 'PP'] <- pi0*bfAll/(pi0*bfAll + (1 - pi0))
    
#    outData <- outData[order(-outData$BF),]
 #   outData$qvalue <- Bayesian.FDR(outData$BF, 1 - pars$pi0)$FDR

    return(outData)
}

##Input: a list of gene sets
###Data

##dnData <-  data[, c("Gene", paste0("dn_", c("damaging", "lof"), "_DD"))]
##mutData <- data[, c("Gene", paste0("mut_", c("damaging", "lof")))]; nFamily = c(4293, 4293)
##geneSetName <- c("fatigue.txt", "rbfox2.txt",  "constrained.txt", "synaptome.txt", "Presynapse.txt", "seizures.txt"); dirGeneSetName = "../../GeneSetAll/List161geneSets/"
#t1 = chooseGeneSet(dnData = dnData, mutData = mutData, nFamily = nFamily,
 #                  geneSetName = geneSetName, dirGeneSetName = dirGeneSetName)

chooseGeneSet <- function(dnData, mutData, nFamily, geneSetName,
                          llkThreshold = 2, dirGeneSetName = "./", outDir = NULL){

    ##Start with no-gene set model (llk0)
    nCol = dim(dnData)[2] - 1
    dD = dnData[, -1]
    muRate = mutData[, -1]
    llk0 <- function(pars){
        hyperGammaMean <- pars[1:nCol]
        hyperBeta <- pars[(nCol + 1):(2*nCol)]
        alpha0 <- tail(pars, 1)
        f0 <- f1 <- 1
        for (j in 1:nCol){
            f0 <- f0*dpois(dD[, j], 2*nFamily[j]*muRate[, j])
            f1<- f1*dnbinom(dD[, j], hyperGammaMean[j]*hyperBeta[j], hyperBeta[j]/(hyperBeta[j] + 2*nFamily[j]*muRate[, j]))
            }
        pi1 <- rep(exp(alpha0)/(1 + exp(alpha0)), dim(dD)[1])
        tllk <- sum(log(f1*pi1 + (1 - pi1)*f0))

        return(-tllk)
    }

    nGeneSet = 0
    pars <- c(2, 30, 1, 1, rep(1, nGeneSet), 1)

    returnValue0 <- nlminb(start = pars, objective=llk0,
       lower = c(rep(1, nCol), rep(0.8, nCol), rep(-10, nGeneSet + 1)),
       upper = c(rep(200, nCol), rep(200, nCol), rep(10, nGeneSet + 1)))

    i0=2
    rValueNew <- NULL
    rValueNew[1] <- returnValue0$objective #value
    ij=1
    indexGene <- 1
    geneOut <- NULL
    tOut <- c(round(returnValue0$objective, 3), "Intercept", "Yes")
    #nMax <- length(p1)


    ##########
    ##Add gene sets
    llk1 <- function(pars){
        hyperGammaMean <- pars[1:nCol]
        hyperBeta <- pars[(nCol + 1):(2*nCol)]
        alpha0 <- pars[(2*nCol+1):(2*nCol + nGeneSet)]
        alphaIntercept <- tail(pars, 1)
        f0 <- f1 <- 1
        for (j in 1:nCol){
            f0 <- f0*dpois(dD[, j], 2*nFamily[j]*muRate[, j])
            f1<- f1*dnbinom(dD[, j], hyperGammaMean[j]*hyperBeta[j], hyperBeta[j]/(hyperBeta[j] + 2*nFamily[j]*muRate[, j]))
            }
        pi1 <- apply(geneSet, 1, function(x){
            #pi1 = exp(geneSet)/(1 + exp(geneSet))
            exp0 <- exp(alphaIntercept + sum(alpha0*x))
            return(exp0/(1 + exp0))})
        tllk <- sum(log(f1*pi1 + (1 - pi1)*f0))

        return(-tllk)
    }

###Add each gene set
    nCol <- ncol(dD)
    nRow <- nrow(dD)
    geneSetM <- matrix(0, ncol = length(geneSetName), nrow = dim(dD)[1])
    for (jj in 1:length(geneSetName)){
        tempGeneSet <- readLines(paste0(dirGeneSetName, geneSetName[jj]))
        geneSetM[dnData[, 1] %in% tempGeneSet, jj] <- 1

    }

    ##Start with the first gene set
    indexGene = 1
    outPar <- NULL
    for (ij in 1:length(geneSetName)){
        geneSet <- data.frame(geneSetM[, sort(unique(c(indexGene, ij)))])
        nGeneSet <- ncol(geneSet)
        message("nGeneSet: ", nGeneSet)
        pars <- c(2, 30, 1, 1, rep(1, nGeneSet), 1)

        returnValue <-    nlminb(start = pars, objective=llk1,
                                 lower = c(rep(1, nCol), rep(0.8, nCol), rep(-10, nGeneSet + 1)),
                                 upper = c(rep(200, nCol), rep(200, nCol), rep(10, nGeneSet + 1)))

        message("logLLK: ", returnValue$objective)
        if ((returnValue$objective - rValueNew[i0-1]) < -llkThreshold){
            indexGene <- c(indexGene, ij)
            rValueNew[i0] <- returnValue$objective #value
            i0 <- i0 + 1
            message("i0: ", i0)
            geneOut <- c(round(returnValue$objective, 2), geneSetName[ij], "Yes")
            tOut <- rbind(tOut, geneOut ) #c(tOut, returnValue$value, pT, returnValue0$par, returnValue$par)
#        tOut <- round(tOut, 3)
            if (!is.null(outDir)){
                write.table(tOut, paste0(outDir, "/TestlogLLKxem.AddGenes.txt"),
                            quote = FALSE, col.names = FALSE, row.names = FALSE)
            }

            parOut <- returnValue$par
            message("ADD gene set: ", geneSetName[ij])

    } else {
        message("REMOVE gene set: ", geneSetName[ij])
        geneOut <- c(round(returnValue$objective, 2), geneSetName[ij], "No")
            tOut <- rbind(tOut, geneOut ) #c(tOut, returnValue$value, pT, returnValue0$par, returnValue$par)
        if (!is.null(outDir)){
                write.table(tOut, paste0(outDir, "/TestlogLLKxem.AddGenes.txt"),
                            quote = FALSE, col.names = FALSE, row.names = FALSE)
            }


    }


    }


    return(list(parGeneSet = parOut, parNoGeneSet = returnValue0$par, geneSetInfo = tOut  ))

}
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
extTADAmcmc <- function(modelName  ,
                         dataDN = NULL, mutRate = NULL, Ndn = NULL,
                  dataCCcase = NULL, dataCCcontrol = NULL, Ncase = NULL, Ncontrol = NULL,
                  geneSet = NULL,
                    nIteration = NULL, nIteration2 = NULL,
                    nThin = NULL, nThin2 = NULL, nCore = 1, nChain = 1,
                         hyperBetaDN0 = NULL,
                         hyperBetaCC0 = NULL,
                         hyper2GammaMeanDN = NULL, hyper2BetaDN = NULL, ##Priors for mean RRs: meanRR ~ gamma(hyper2GammaMeanDN, hyper2BetaDN)
                         hyper2GammaMeanCC = NULL, hyper2BetaCC = NULL,
                  alpha01Mean = 0, alpha01SD = 2, alpha01U = 2,

                         upperPi0 = 0.5, lowerPi0 = 0, lowerBeta = 0, ##Lower and upper limits of pi: should be default
                         lowerHyperGamma = 1, lowerGamma = 1, #Should be default
                         betaPars = c(6.7771073, -1.7950864, -0.2168248), #Adjust beta's values: should be default
                    adjustHyperBeta = as.integer(1), ##1 if want to adjust beta, 0 if don't want to adjust beta
                  autoAdjustHyperBeta = FALSE,
                  sigmaPrior = 2)
     {

         if (is.null(nIteration))
             stop("======\nNo input for the nIteration parameter\n=====")
         if ((is.null(dataDN) | is.null(mutRate)) & (is.null(dataCCcase) | is.null(dataCCcontrol)))
             stop("Need to have input data: only DN, only CC or DN + CC")
         if (is.null(nThin))
             nThin = ifelse(nIteration > 1000, floor(nIteration/1000), 1)
         if (nCore != nChain)
             warning("nCore is different from nChain")
         if (!is.null(dataDN))
             Ngene <- dim(dataDN)[1]
         if (!is.null(dataCCcase))
             Ngene <- dim(dataCCcase)[1]
         message("\nThere are ", Ngene, " genes in this analysis")

         if (is.null(hyper2GammaMeanDN) & !is.null(dataDN))
             hyper2GammaMeanDN <- rep(1, dim(dataDN)[2])
         if (is.null(hyper2BetaDN) & !is.null(dataDN))
             hyper2BetaDN <- rep(0.025, dim(dataDN)[2])
        if (is.null(hyper2GammaMeanCC) & !is.null(dataCCcase))
             hyper2GammaMeanCC <- rep(1, dim(dataCCcase)[2])
         if (is.null(hyper2BetaCC) & !is.null(dataCCcase))
             hyper2BetaCC <- rep(0.2, dim(dataCCcase)[2])
         if ((adjustHyperBeta == 0) & is.null( hyperBetaDN0))
             hyperBetaDN0 <- rep(1, dim(dataDN)[2])
         if ((adjustHyperBeta == 0) & is.null( hyperBetaCC0) & !is.null(dataCCcase))
             hyperBetaCC0 <- rep(1, dim(dataCCcase)[2])

         NCdn <- dim(dataDN)[2]
         NCcc <- dim(dataCCcase)[2]
         message("\nNCdn: ", NCdn, "; NCcc: ", NCcc)
         
         if (is.null(Ncase))              Ncase = 0
         if (is.null(Ncontrol)) Ncontrol = 1
         if (is.null(Ndn)) Ndn = 0
         if (is.null(hyperBetaDN0)) hyperBetaDN0 <- rep(1, length(dataDN[1, ]))
         if (is.null(hyperBetaCC0)) hyperBetaCC0 <- rep(4, length(dataCCcase[1, ]))
         if (is.null(hyper2GammaMeanDN)) hyper2GammaMeanDN = rep(1, length(dataDN[1, ]))
         if (is.null(hyper2GammaMeanCC)) hyper2GammaMeanCC = rep(4, length(dataCCcase[1, ]))
         if (is.null(hyper2BetaDN)) hyper2BetaDN = rep(0.05, length(dataDN[1, ]))
         if (is.null(hyper2BetaCC)) hyper2BetaCC = rep(0.2, length(dataCCcase[1, ]))


         Tgs <- 0
         if (!is.null(geneSet)) {
             Tgs = dim(geneSet)[2]
             Ngs = Tgs
         } else {
             geneSet <- data.frame(sample(0:1, Ngene, replace = TRUE))
             Ngs = dim(geneSet)[2]
             message("No gene set input\n")

         }


         modelData <- list(NN = Ngene, #Gene numbers
                           K = 2, #Hypothesis numbers: should be default
                           NCdn = NCdn, #Number of de novo classes
                           NCcc = NCcc,
                           Ndn = array(Ndn), # Family numbers
                           Ncase = array(Ncase), Ncontrol = array(Ncontrol), Ntotal = array(Ncase + Ncontrol),
                           dataDN = dataDN, # Denovo data
                           mutRate = mutRate, # Mutation rates
                           dataCCcase = data.frame(dataCCcase), dataCCtotal = data.frame(dataCCcase + dataCCcontrol),
                           thetaH0 = array(Ncase/(Ncase + Ncontrol)),
                           betaPars = array(betaPars), #Adjust beta's values: should be default
                           adjustHyperBeta = adjustHyperBeta, ##1 if want to adjust beta, 0 if don't want to adjust beta
                           upperPi0 = upperPi0, lowerPi0 = lowerPi0, lowerBeta = lowerBeta, ##Lower and upper limits of pi: should be default
                          lowerHyperGamma = lowerHyperGamma, lowerGamma = lowerGamma, #Should be default
                          hyperBetaDN0 = array(hyperBetaDN0),
                          hyperBetaCC0 = array(hyperBetaCC0),
                            hyper2GammaMeanDN = array(hyper2GammaMeanDN),##Priors for mean RRs: meanRR ~ gamma(hyper2GammaMeanDN, hyper2BetaDN)
                          hyper2BetaDN = array(hyper2BetaDN),
                          hyper2GammaMeanCC = array(hyper2GammaMeanCC),
                          hyper2BetaCC = array(hyper2BetaCC),
                          Tgs = Tgs,
                          Ngs = Ngs,
                          geneSet = geneSet,
                          sigmaPrior = sigmaPrior,
                          alpha01Mean = alpha01Mean, alpha01SD = alpha01SD,
                          alpha01U = alpha01U)

         message("\n=============FIRST TIME ===============\n")
         message("\nSampling with nter = ", nIteration, " and nThin = ", nThin, "\n")
         message("\nThe model ", deparse(substitute(modelName)), " is used\n")
         message("\nNumber of gene sets: ", dim(modelData$geneSet)[2], "\n")
         mcmcModel <- stan(model_code = modelName,
                        data = modelData, ##Model data as described above
                        iter = nIteration, chains = nChain, cores = nCore, thin = nThin)

         ####################Finish the FIRST TIME
####Re-sample using new hyper betas
         if ( autoAdjustHyperBeta){
             if (!is.null(nIteration2))
                 nIteration = nIteration2

             if (is.null(nThin2))
                 nThin2 = ifelse(nIteration > 1000, floor(nIteration/1000), 1)
             nThin = nThin2
             mcmcData <- as.data.frame(mcmcModel)
             cName <- colnames(mcmcData)
         hyperGammaMeanNameCC <- cName[grep("hyperGammaMeanCC", cName)]
         hyperGammaMeanNameDN <- cName[grep("hyperGammaMeanDN", cName)]

             hyperGammaMeanDN0 <- as.numeric(apply(data.frame(mcmcData[, hyperGammaMeanNameDN]), 2, median))
             hyperGammaMeanCC0 <- as.numeric(apply(data.frame(mcmcData[, hyperGammaMeanNameCC]), 2, median))

             hyperBetaDN0 <- exp(betaPars[1]*hyperGammaMeanDN0^(betaPars[2]) + betaPars[3])
             hyperBetaCC0 <- exp(betaPars[1]*hyperGammaMeanCC0^(betaPars[2]) + betaPars[3])

             adjustHyperBeta <- 0

         modelData <- list(NN = Ngene, #Gene numbers
                           K = 2, #Hypothesis numbers: should be default
                           NCdn = NCdn, #Number of de novo classes
                           NCcc = NCcc,
                           Ndn = array(Ndn), # Family numbers
                           Ncase = array(Ncase), Ncontrol = array(Ncontrol), Ntotal = array(Ncase + Ncontrol),
                           dataDN = dataDN, # Denovo data
                           mutRate = mutRate, # Mutation rates
                           dataCCcase = data.frame(dataCCcase), dataCCtotal = data.frame(dataCCcase + dataCCcontrol),
                           thetaH0 = array(Ncase/(Ncase + Ncontrol)),
                           betaPars = array(betaPars), #Adjust beta's values: should be default
                           adjustHyperBeta = adjustHyperBeta, ##1 if want to adjust beta, 0 if don't want to adjust beta
                           upperPi0 = upperPi0, lowerPi0 = lowerPi0, lowerBeta = lowerBeta, ##Lower and upper limits of pi: should be default
                          lowerHyperGamma = lowerHyperGamma, lowerGamma = lowerGamma, #Should be default
                          hyperBetaDN0 = array(hyperBetaDN0),
                           hyperBetaCC0 = array(hyperBetaCC0),
                            hyper2GammaMeanDN = array(hyper2GammaMeanDN),##Priors for mean RRs: meanRR ~ gamma(hyper2GammaMeanDN, hyper2BetaDN)
                          hyper2BetaDN = array(hyper2BetaDN),
                          hyper2GammaMeanCC = array(hyper2GammaMeanCC),
                          hyper2BetaCC = array(hyper2BetaCC),
                          alpha01Mean = alpha01Mean, alpha01SD = alpha01SD)

             message("\n=============SECOND TIME ===============\n")
             message("\nThis process is using beta values estimated from the previous process\n")
         message("\nSampling with nter = ", nIteration, " and nThin = ", nThin, "\n")
         message("\nThe model ", deparse(substitute(modelName)), " is used\n")
         mcmcModel <- stan(model_code = modelName,
                        data = modelData, ##Model data as described above
                        iter = nIteration, chains = nChain, cores = nCore, thin = nThin)
         }



         return(mcmcModel)
          }



###Bayes Factor For Case-control

##########There are three models written in this script: DNandCCextTADA (de novo + case control), DNextTADA (only de novo), CCextTADA (only case-control)
##########Users can change prior information as they wish
library("rstan")
library("coda")
########################################
#######De novo + Case control
#June 1, 2017: Hoang Nguyen

DNandCCextTADA <- "
data {
    int<lower=1> NN; //Number of genes
    int<lower=1> K; //Number of classes
    int<lower=1> NCdn; //Number of de novo classes
    int<lower=1> NCcc; //Number of case-control classes
    int Ndn[NCdn]; //Number of trios
    int Ncase[NCcc]; //Number of cases
    int Ncontrol[NCcc]; //Number of controls
    int Ntotal[NCcc]; //Number of cases + controls

    int dataDN[NN, NCdn]; //denovo data: Kdn classes

    real mutRate[NN, NCdn]; //mutation rates: Kdn classes
    int dataCCcase[NN, NCcc]; //case-control data: Kcc classes
    int dataCCtotal[NN, NCcc]; //case+control data: Kcc classes
    real betaPars[3]; //These parameters are used to adjust beta values as a non-linear function of gamma means

    real<lower=0> thetaH0[NCcc]; // Ncase/(Ncase + Ncontrol)

    //These below parameters should be default
//    real<lower=0> upperPi0;
//    real<lower=0> lowerPi0;
    real<lower=0> lowerHyperGamma; //Low limit for mean of relative risks
    real<lower=0> lowerGamma; //Low limit for relative risks
    real<lower=0> lowerBeta;
    real<lower=0> hyperBetaDN0[NCdn];
    real<lower=0> hyperBetaCC0[NCcc];
    int<lower=0> adjustHyperBeta; //Option to adjust betas or not; if this option is used (!=0) then beta is a non-linear function of gamma means

    int<lower=0> Tgs; ##Test gene sets
    int<lower=0> Ngs;
    real geneSet[NN, Ngs];

    real alpha01Mean;
    real<lower=0> alpha01SD;
    real alpha01U;

    }

parameters {
//    real<lower=lowerPi0,upper=upperPi0> pi0; //Proportion of risk genes
    real<upper=alpha01U> alpha01;
    real alpha0[Tgs];
    real<lower=lowerHyperGamma> hyperGammaMeanDN[NCdn]; //Hyper parameters for de novo relative risks
    real<lower=lowerGamma> gammaMeanDN[NCdn]; //parameters (in the sampling process) for de novo relative risks

//    real<lower=lowerBeta> betaDN[Kdn]; //Rate parameters for de novo relative risks
    //real<lower=lowerBeta> hyperBetaDN[Kdn]; //Hyper rate parameters for de novo relative risks
    real<lower=lowerHyperGamma> hyperGammaMeanCC[NCcc]; //Hyper parameters for case-control relative risks
    //ordered[Kcc] hyperGammaMeanCC;
    real<lower=lowerGamma> gammaMeanCC[NCcc]; //parameters (in the sampling process) for de novo relative risks: gammaMeanCC ~ gamma(hyperGammaMeanCC*hyperBetaCC, hyperBetaCC)


}

transformed parameters {

    real sumGeneSet[NN];
    real pi0[NN]; //Proportion of risk genes
    real hyperBetaCC[NCcc];
    real hyperBetaDN[NCdn];

    for (i in 1:NN){
      sumGeneSet[i] = alpha01;
      if (Tgs > 0){
        for (j in 1:Ngs){
          sumGeneSet[i] = sumGeneSet[i] + alpha0[j]*geneSet[i, j] ;
        }}

    }
    for (i in 1:NN){
           pi0[i] = exp(sumGeneSet[i])/(1 + exp(sumGeneSet[i]));
     }



    if (adjustHyperBeta != 0) {
      for (i1i in 1:NCcc){
            hyperBetaCC[i1i] = exp(betaPars[1]*hyperGammaMeanCC[i1i]^(betaPars[2]) + betaPars[3]); //,  hyperBetaMax);

       }

      for (i2i in 1:NCdn){
            hyperBetaDN[i2i] = exp(betaPars[1]*hyperGammaMeanDN[i2i]^(betaPars[2]) + betaPars[3]); //,  hyperBetaMax);

       }
   }
    else {
        hyperBetaCC = hyperBetaCC0;
        hyperBetaDN = hyperBetaDN0;
        }
    }
 model {
     int newIndex;
     real ps[K];
    alpha01 ~ normal(alpha01Mean, alpha01SD);
    if (Tgs > 0)
       for (k in 1:Ngs)
         alpha0[k] ~ normal(0, 2);


     //Both CC + DN

     //Case control data: sample for hyper priors (NPcc populations and Kcc categories)
     for (ip in 1:NCcc){
         hyperGammaMeanCC[ip]  ~ gamma(1, 0.1); //gamma(1, 0.05); //normal(15, 10); //gamma(1, 0.1);

     }

     for (ip in 1:NCcc){
         gammaMeanCC[ip] ~ gamma(hyperGammaMeanCC[ip]*hyperBetaCC[ip], hyperBetaCC[ip]);
         }

  //De novo data: sample for hyper priors (NPdn populations and Kdn categories)
    for (ip in 1:NCdn){
          hyperGammaMeanDN[ip]	~ gamma(1, 0.02); //gamma(1, 0.05); //gamma(1, 0.1); //normal(15, 10);
    }

    for (ip in 1:NCdn){
          gammaMeanDN[ip] ~ gamma(hyperGammaMeanDN[ip]*hyperBetaDN[ip], hyperBetaDN[ip]);
           }

////Main program
//Loop through data points
////
     for (ii in 1:NN){
         ps[1] = log1m(pi0[ii]);
         ps[2] = log(pi0[ii]);
   //For de novo data
         for (jj in 1:NCdn){
             ps[1] = ps[1] + poisson_lpmf(dataDN[ii, jj] | Ndn[jj]*2*mutRate[ii, jj]); //Add Null hypothesis
             ps[2] = ps[2] + poisson_lpmf(dataDN[ii, jj] | Ndn[jj]*2*mutRate[ii, jj]*gammaMeanDN[jj]); //Add alternative hypothesis
             }
    //For case-control data
         for (jj in 1:NCcc){
             ps[1] = ps[1] + binomial_lpmf(dataCCcase[ii, jj] | dataCCtotal[ii, jj], thetaH0[jj]); //Add Null hypothesis
             ps[2] = ps[2] + binomial_lpmf(dataCCcase[ii, jj] | dataCCtotal[ii, jj], gammaMeanCC[jj]*Ncase[jj]/(gammaMeanCC[jj]*Ncase[jj] + Ncontrol[jj])); //Add alternative hypothesis
             }
         target += log_sum_exp(ps);
         //increment_log_prob(log_sum_exp(ps));
         }
}
"


########################################
#######De novo only
DNextTADA <- "
data {
    int<lower=1> NN; //Number of genes
    int<lower=1> K; //Number of classes
    int<lower=1> NCdn; //Number of de novo classes
    int Ndn[NCdn]; //Number of trios

    int dataDN[NN, NCdn]; //denovo data: Kdn classes
    real mutRate[NN, NCdn]; //mutation rates: Kdn classes
    real betaPars[3]; //These parameters are used to adjust beta values as a non-linear function of gamma means

    //These below parameters should be default
    real<lower=0> upperPi0;
    real<lower=0> lowerPi0;
    real<lower=0> lowerHyperGamma; //Low limit for mean of relative risks
    real<lower=0> lowerGamma; //Low limit for relative risks
    real<lower=0> lowerBeta;
    real<lower=0> hyperBetaDN0[NCdn];
    int<lower=0> adjustHyperBeta; //Option to adjust betas or not; if this option is used (!=0) then beta is a non-linear function of gamma means
    int<lower=0> Tgs; //Gene-set numbers (Tgs = 0 if no gene set is used)
    int<lower=0> Ngs; //Gene-set numbers: should be from input gene sets
    real geneSet[NN, Ngs];
    real<lower=0> sigmaPrior;
    real alpha01Mean;
    real<lower=0> alpha01SD;
    real alpha01U;

    }

parameters {
    real<upper=alpha01U> alpha01;
    real alpha0[Tgs];
    real<lower=lowerHyperGamma> hyperGammaMeanDN[NCdn]; //Hyper parameters for de novo relative risks
    real<lower=lowerGamma> gammaMeanDN[NCdn]; //parameters (in the sampling process) for de novo relative risks

}

transformed parameters {
    real pi0[NN]; //Proportion of risk genes
    real sumGeneSet[NN];
    real hyperBetaDN[NCdn];

    for (i in 1:NN){
      sumGeneSet[i] = alpha01;
      if (Tgs > 0){
         for (j in 1:Ngs){
         sumGeneSet[i] = sumGeneSet[i] + alpha0[j]*geneSet[i, j] ;
       }}
    }
    for (i in 1:NN)
       pi0[i] =   exp(sumGeneSet[i])/(1 + exp(sumGeneSet[i]));

    if (adjustHyperBeta != 0) {
      for (i2i in 1:NCdn){
            hyperBetaDN[i2i] = exp(betaPars[1]*hyperGammaMeanDN[i2i]^(betaPars[2]) + betaPars[3]);

       }
   }
    else {
        hyperBetaDN = hyperBetaDN0;
        }
    }
 model {
     int newIndex;
     real ps[K];
     real sumDN[2];
//     pi0 ~ beta(1, 5); //prior for the proportion of risk genes
    alpha01 ~ normal(alpha01Mean, alpha01SD);
    if (Tgs > 0)
      for (k in 1:Ngs)
        alpha0[k] ~ normal(0, sigmaPrior); //normal(0, 2);
     //Both CC + DN


  //De novo data: sample for hyper priors (NPdn populations and Kdn categories)
    for (ip in 1:NCdn){
          hyperGammaMeanDN[ip]	~ gamma(1, 0.05); //gamma(1, 0.1); //normal(15, 10);
    }

    for (ip in 1:NCdn){
          gammaMeanDN[ip] ~ gamma(hyperGammaMeanDN[ip]*hyperBetaDN[ip], hyperBetaDN[ip]);
           }

////Main program
//Loop through data points
////
     for (ii in 1:NN){
         ps[1] = log1m(pi0[ii]);
         ps[2] = log(pi0[ii]);
   //For de novo data
         for (jj in 1:NCdn){
             ps[1] = ps[1] + poisson_lpmf(dataDN[ii, jj] | Ndn[jj]*2*mutRate[ii, jj]); //Add Null hypothesis
             ps[2] = ps[2] + poisson_lpmf(dataDN[ii, jj] | Ndn[jj]*2*mutRate[ii, jj]*gammaMeanDN[jj]); //Add alternative hypothesis
             }
         target += log_sum_exp(ps);
         //increment_log_prob(log_sum_exp(ps));
         }
}
"


################CC model
########################

CCextTADA <- "
data {
    int<lower=1> NN; //Number of genes
    int<lower=1> K; //Number of classes
    int<lower=1> NCcc; //Number of case-control classes
    int Ncase[NCcc]; //Number of cases
    int Ncontrol[NCcc]; //Number of controls
    int Ntotal[NCcc]; //Number of cases + controls

    int dataCCcase[NN, NCcc]; //case-control data: Kcc classes
    int dataCCtotal[NN, NCcc]; //case+control data: Kcc classes
    real betaPars[3]; //These parameters are used to adjust beta values as a non-linear function of gamma means

    real<lower=0> thetaH0[NCcc]; // Ncase/(Ncase + Ncontrol)

    //These below parameters should be default
//    real<lower=0> upperPi0;
//    real<lower=0> lowerPi0;
    real<lower=0> lowerHyperGamma; //Low limit for mean of relative risks
    real<lower=0> lowerGamma; //Low limit for relative risks
    real<lower=0> lowerBeta;
    real<lower=0> hyperBetaCC0[NCcc];
    int<lower=0> adjustHyperBeta; //Option to adjust betas or not; if this option is used (!=0) then beta is a non-linear function of gamma means

    int<lower=0> Tgs; ##Test gene sets
    int<lower=0> Ngs;
    real geneSet[NN, Ngs];
    real alpha01Mean;
    real<lower=0> alpha01SD;
    real alpha01U;

    }

parameters {
//    real<lower=lowerPi0,upper=upperPi0> pi0; //Proportion of risk genes
    real<upper=alpha01U> alpha01; 
    real alpha0[Tgs];
    real<lower=lowerHyperGamma> hyperGammaMeanCC[NCcc]; //Hyper parameters for case-control relative risks
    //ordered[Kcc] hyperGammaMeanCC;
    real<lower=lowerGamma> gammaMeanCC[NCcc]; //parameters (in the sampling process) for de novo relative risks: gammaMeanCC ~ gamma(hyperGammaMeanCC*hyperBetaCC, hyperBetaCC)


}

transformed parameters {

    real sumGeneSet[NN];
    real pi0[NN]; //Proportion of risk genes
    real hyperBetaCC[NCcc];


    for (i in 1:NN){
      sumGeneSet[i] = alpha01;
      if (Tgs > 0){
        for (j in 1:Ngs){
          sumGeneSet[i] = sumGeneSet[i] + alpha0[j]*geneSet[i, j] ;
        }}

    }
    for (i in 1:NN){
           pi0[i] = exp(sumGeneSet[i])/(1 + exp(sumGeneSet[i]));
     }



    if (adjustHyperBeta != 0) {
      for (i1i in 1:NCcc){
            hyperBetaCC[i1i] = exp(betaPars[1]*hyperGammaMeanCC[i1i]^(betaPars[2]) + betaPars[3]); //,  hyperBetaMax);

       }

   }
    else {
        hyperBetaCC = hyperBetaCC0;
        }
    }
 model {
     int newIndex;
     real ps[K];
    alpha01 ~ normal(alpha01Mean, alpha01SD);
    if (Tgs > 0)
       for (k in 1:Ngs)
         alpha0[k] ~ normal(0, 2);


     //Case control data: sample for hyper priors (NPcc populations and Kcc categories)
     for (ip in 1:NCcc){
         hyperGammaMeanCC[ip]  ~ gamma(1, 0.1); //gamma(1, 0.05); //normal(15, 10); //gamma(1, 0.1);

     }

     for (ip in 1:NCcc){
         gammaMeanCC[ip] ~ gamma(hyperGammaMeanCC[ip]*hyperBetaCC[ip], hyperBetaCC[ip]);
         }


////Main program
//Loop through data points
////
     for (ii in 1:NN){
         ps[1] = log1m(pi0[ii]);
         ps[2] = log(pi0[ii]);
    //For case-control data
         for (jj in 1:NCcc){
             ps[1] = ps[1] + binomial_lpmf(dataCCcase[ii, jj] | dataCCtotal[ii, jj], thetaH0[jj]); //Add Null hypothesis
             ps[2] = ps[2] + binomial_lpmf(dataCCcase[ii, jj] | dataCCtotal[ii, jj], gammaMeanCC[jj]*Ncase[jj]/(gammaMeanCC[jj]*Ncase[jj] + Ncontrol[jj])); //Add alternative hypothesis
             }
         target += log_sum_exp(ps);
         //increment_log_prob(log_sum_exp(ps));
         }
}
"


extTADAsimulator <- function(N, muAll, 
                             pi0,
                             genesetPar = NULL, ## Gene set parameters
                       beta0 = NULL, ## Mutation-rate parameter
                       gamma.mean.dn, beta.dn,
                       gamma.mean.CC, beta.CC,
                       rho1, nu1, rho0, nu0, tradeoff=FALSE,
                       adjustMutationRate = FALSE) {
      m <- dim(muAll)[1] # number of genes
      K <- dim(muAll)[2] # number of mutational categories
      pp0 <- pi0
      z <- sample(c(1, 0), size = m, replace = TRUE, prob = c(pp0, 1 - pp0))
        gamma.dn <- array(1, dim=c(m,K))
        gamma.CC <- array(1, dim=c(m,K))
        q <- array(0, dim=c(m,K))
        x <- array(0, dim=c(m,3*K))
        k <- sum(z==1)

      if (!is.null(genesetPar)){
          geneSet <- rep(0, length(z))
          geneSet[z == 1] <- sample(c(1, 0), length(z[z == 1]), replace = TRUE, prob = c(genesetPar[1], 1 - genesetPar[1]))
          geneSet[z == 1] <- sample(c(1, 0), length(z[z == 1]), replace = TRUE, prob = c(genesetPar[2], 1 - genesetPar[2]))          
      }
      
      for (j in 1:K) {
                mu <- muAll[, j]
                    # sample de novo
                    gamma.dn[z==1, j] <- rgamma(k, gamma.mean.dn[j]*beta.dn[j], beta.dn[j])
                    col <- 3*(j-1)+1
                    x[,col] <- rpois(m, 2 * mu  * gamma.dn[,j] * N$dn[j])

                    # sample case-control
                    gamma.CC[z==1, j] <- rgamma(k, gamma.mean.CC[j]*beta.CC[j], beta.CC[j])
                    q[z==0, j] <- rgamma(m-k, rho0[j], nu0[j])
                    if (tradeoff==FALSE) {
                              q[z==1, j] <- rgamma(k, rho1[j], nu1[j])
                    } else {
                              q[z==1, j] <- mu[z==1] / (delta[j] * gamma.CC[z==1, j])
                    }
                    x[,col+1] <- rpois(m, q[,j] * gamma.CC[,j] *N$ca[j])
                    x[,col+2] <- rpois(m, q[,j] * N$cn[j])

        }

        sample.info <- cbind(z, gamma.dn, gamma.CC, q, x, muAll, geneSet)

        return (list(sample=x, sample.info=sample.info))
}
gTADAfromBF <- "
data {
    int<lower=1> NN; //Number of genes

    //These below parameters should be default
    int<lower=0> Ngs; //Gene-set numbers: should be from input gene sets
    real geneSet[NN, Ngs];
    real<lower=0> sigmaPrior;
    real alpha01Mean;
    real<lower=0> alpha01SD;
    real<lower=0> bf0[NN];
    real<lower=0> pH0[NN];

    }

parameters {
    real alpha0[Ngs + 1];

}

transformed parameters {
    real pi0[NN]; //Proportion of risk genes
    real sumGeneSet[NN];

    for (i in 1:NN){
      sumGeneSet[i] = alpha0[1];
      for (j in 1:Ngs){
         sumGeneSet[i] = sumGeneSet[i] + alpha0[j+1]*geneSet[i, j] ;
       }
    }

    for (i in 1:NN)
       pi0[i] =   exp(sumGeneSet[i])/(1 + exp(sumGeneSet[i]));
}
 model {
//     int newIndex;
     real ps[2];
//     pi0 ~ beta(1, 5); //prior for the proportion of risk genes
    alpha0[1] ~ normal(alpha01Mean, alpha01SD);
    for (k in 1:Ngs)
        alpha0[k+1] ~ normal(0, sigmaPrior); //normal(0, 2);
    
////Main program
//Loop through data points
////
     for (ii in 1:NN){
         ps[1] = log(pi0[ii]) + log(bf0[ii]) + log(pH0[ii]); //H1: risk-gene hypothesis
         ps[2] = log1m(pi0[ii]) + log(pH0[ii]); //H0: non-risk gene hypothesis
         target += log_sum_exp(ps);
         //increment_log_prob(log_sum_exp(ps));
         }
}
"

pH0.model.CC <- function(data.cc, n.cc, rho0 = NULL, nu0 = NULL){
    p.null.CC <- matrix(0, nrow = dim(data.cc)[1], ncol = dim(data.cc)[2]/2)
    for (j in 1:floor(dim(data.cc)[2]/2)){
        print(head(data.cc))
        dataTemp <- data.cc[, paste0(c("cc_case", "cc_control"), j), drop = FALSE]
        if (is.null(nu0)){
            e.nu <- 200
            } else e.nu = nu0[j]
        if (is.null(rho0)){
            e.rho <- e.nu*mean(rowSums(dataTemp))/(n.cc$ncase[j] + n.cc$ncontrol[j])
            } else e.rho = rho0[j]
###NULL model for CC data
        pTemp0 <- apply(dataTemp, 1, function(y){
            evidence.null.CC(x = list(ca = y[1], cn = y[2]), N = list(ca = n.cc$ncase[j], cn = n.cc$ncontrol[j]),
                             rho0 = e.rho, nu0 = e.nu)$total})
        p.null.CC[, j] <- pTemp0
        
    }
    return(p.null.CC)
    
}
#######De novo data
pH0.model.DN <- function(data.dn, data.mut, n.dn){
    p.null.DN <- matrix(0, nrow = dim(data.dn)[1], ncol = dim(data.dn)[2])
    ##Checking the names of DNMs and mutation rates: they should be the same
    dn.names <- sapply(strsplit(colnames(data.dn), split = "_", fixed = TRUE), function(x) x[2])
    mut.names <- sapply(strsplit(colnames(data.mut), split = "_", fixed = TRUE), function(x) x[2])
    if (!identical(dn.names, mut.names))
        warning('\nNames of de novo mutations and mutations are not identical\nThis can result in ')
    
    for (j in 1:floor(dim(data.dn)[2])){
        p.null.DN[, j] <- dpois(x = data.dn[, j], lambda = 2*n.dn[j]*data.mut[, j])
        }
    return(p.null.DN)
    
}


################
###Add gene set
llk1 <- function(pars){
        hyperGammaMean <- pars[1:nCol]
        hyperBeta <- pars[(nCol + 1):(2*nCol)]
        alpha0 <- pars[(2*nCol+1):(2*nCol + nGeneSet)]
        alphaIntercept <- tail(pars, 1)
        f0 <- f1 <- 1
        for (j in 1:nCol){
            f0 <- f0*dpois(dD[, j], 2*nFamily[j]*muRate[, j])
            f1<- f1*dnbinom(dD[, j], hyperGammaMean[j]*hyperBeta[j], hyperBeta[j]/(hyperBeta[j] + 2*nFamily[j]*muRate[, j]))
            }
        pi1 <- apply(geneSet, 1, function(x){
            #pi1 = exp(geneSet)/(1 + exp(geneSet))
            exp0 <- exp(alphaIntercept + sum(alpha0*x))
            return(exp0/(1 + exp0))})
        tllk <- sum(log(f1*pi1 + (1 - pi1)*f0))

        return(-tllk)
    }

llk0 <- function(pars){
        hyperGammaMean <- pars[1:nCol]
        hyperBeta <- pars[(nCol + 1):(2*nCol)]
        alpha0 <- pars[(2*nCol+1):(2*nCol + nGeneSet + 1)]
#        alphaIntercept <- tail(pars, 1)

        f0 <- f1 <- 1
        for (j in 1:nCol){
            f0 <- f0*dpois(dD[, j], 2*nFamily[j]*muRate[, j])
            f1<- f1*dnbinom(dD[, j], hyperGammaMean[j]*hyperBeta[j], hyperBeta[j]/(hyperBeta[j] + 2*nFamily[j]*muRate[, j]))
            }
        pi1 <- rep(exp(alpha0)/(1 + exp(alpha0)), dim(dD)[1])

        tllk <- sum(log(f1*pi1 + (1 - pi1)*f0))

        return(-tllk)
    }

##The main function is here; some useful function can be found below the 'plotParHeatmap' function
plotParHeatmap <- function(pars, mcmcResult, nThin = NULL, color = "blue",
                           xLim = NULL, yLim = NULL,
                           mainLab = NULL, xLab = NULL, yLab = NULL,
                           maxk0 = 500, cprob = c(0.5, 0.05)){
    mcmcDataFrame <- as.data.frame(mcmcResult)
    if (!is.null(nThin))
        mcmcDataFrame <- mcmcDataFrame[seq(1, dim(mcmcDataFrame)[1], by = nThin),]
    allPars <- colnames(mcmcDataFrame)

    x <- mcmcDataFrame[, pars[1]]
    y <- mcmcDataFrame[, pars[2]]
    sc1<-sd(x)
    sc2<-sd(y)
    fit <- locfit(~x+y,scale=c(sc1,sc2), maxk = maxk0)
    lev <- sort(fitted(fit))[floor(cprob*length(x))]
    max.i <- fitted(fit)==max(fitted(fit))
    mode.x<-x[max.i][[1]]
    mode.y<-y[max.i][[1]]
    if (is.null(yLim))
        yLim <- c(0, max(y))
    if (is.null(xLim))
        xLim <- c(min(x), max(x))
    if (is.null(xLab))
        xLab <- pars[1]
    if (is.null(yLab))
        yLab <- pars[2]

    plot(mode.x,mode.y,
     #xaxt="n",
         main = mainLab, #diseaseName,
         xlim = xLim, ylim = yLim,
         pch=-1, xlab = xLab, ylab = yLab)

my.color.palette=colorRampPalette(c("white", color), space = "Lab")
plot(fit, type="image", m=300, add=TRUE,
     col=my.color.palette(25)[-c(2,3,5,6,8,10,12)])

plot(fit,add=TRUE,lev=lev,m=100,label=paste(as.character(100*(1-cprob)),"%",sep=""), labcex=0.75,vfont=c("sans serif","bold"))
points(mode.x,mode.y,pch=3)

}



##################
################


library(locfit)
loc2plot <- function(x,y,cprob=0.5, xlim,gxlim,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,scale=c(sc1,sc2), maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	else fit <- locfit(~x+y,scale=c(sc1,sc2), xlim=xlim,maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	lev <- sort(fitted(fit))[floor(cprob*length(x))]
	plot(fit,lev=lev,m=100,label=paste(as.character(100*(1-cprob)),"%",sep=""),
	xlim=gxlim,...)
}

# finds the mode for a bivariate density
loc2mode <- function(x,y,alpha=0.5,xlim,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2), maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),xlim=xlim,maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	tt <- max(fitted(fit))
	wt <- fitted(fit) == tt
	c(x[wt][1],y[wt][1])
}

# this works for univariate data; gives a vector with
# mode (global), hpd1_low, hpd1_high, hpd2_low, hpd2_hi, etc
#The reason for multiple hpd limits is if you have multimodal data.
#prob = prob *outside* the limit; i.e for a normal distribution 0.05 is expected to give
#     the 0.025 and 0.975 quantiles.
# this won't work for weighted data, use loc1statsx instead.
# xlim is optional - use it to define the support of the density.
loc1stats <- function(x,xlim,prob=0.05,...)
{
	if(missing(xlim)){
		fit <- locfit(~x)
	}
	else {
		fit <- locfit(~x,xlim=xlim)
	}
	fx <- fitted(fit)
	x.modef <- max(fx)
	x.mode <- x[fx == x.modef]
	if(!missing(xlim)){
		if(predict(fit,xlim[1]) > x.modef){
			x.modef <- predict(fit,xlim[1])
			x.mode <- xlim[1]
		}
		if(predict(fit,xlim[2]) > x.modef){
			x.modef <- predict(fit,xlim[2])
			x.mode <- xlim[2]
		}
	}

	if(length(x.mode)>1)x.mode <- x.mode[1]
	lev <- sort(fx)[floor(prob*length(x))]
#	print("log difference from max is ")
#	print(log(x.modef)-log(lev))
	l1 <- list()
	l1[[1]] <- x.mode
	indx <- order(x)
	ii <- 2
	flip <- TRUE
	for(j in 2:length(x)){
		if(flip && fx[indx[j]] > lev){
			l1[[ii]] <- x[indx[j-1]]
			if(j==2 && !missing(xlim)){
				if(predict(fit,xlim[1]) >= lev)l1[[ii]] <- xlim[1]
			}
			flip <- FALSE
			ii <- ii+1
		}
		else if(!flip && fx[indx[j]] < lev){
			l1[[ii]] <- x[indx[j]]
			flip <- TRUE
			ii <- ii+1
		}
		if(!flip && !missing(xlim) && j == length(x)){
			l1[[ii]] <- xlim[2]
			flip <- TRUE
		}
	}
	if(!flip)stop("HPD interval not closed")
	as.numeric(l1)
}


