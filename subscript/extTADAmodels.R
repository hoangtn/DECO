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


