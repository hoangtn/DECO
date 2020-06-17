-   [Introduction](#introduction)
-   [Usage](#usage)
    -   [Run `DECO` for both variant and gene-set
        data.](#run-textttdeco-for-both-variant-and-gene-set-data.)
        -   [Data format](#data-format)
        -   [Use `DECO` function](#run1)
        -   [Outputs of `DECO`](#outputs-of-textttdeco)
    -   [Using genetic parameters from previous studies or from 1) to
        run
        `DECO`](#using-genetic-parameters-from-previous-studies-or-from-1-to-run-textttdeco)
        -   [Read a gene-set (e.g. pLI09 in this
            example)](#read-a-gene-set-e.g.pli09-in-this-example)
        -   [Run `DECO`](#run-textttdeco)
        -   [Results](#results)
-   [Citation](#citation)

Introduction
============

`DECO` jointly analyzes de novo mutations (DNMs) of parent-offspring
trios, rare case/control variants and gene-set information to:

1.  conduct gene-set enrichment analysis (GSEA),

2.  prioritize risk genes for the tested disease.

Usage
=====

Users can go directly to [this step](#run1) to run `DECO`. However, some
details are described here.

There are two ways to obtain analysis results of `DECO`:

1.  Running `DECO` for both variant and gene-set data. It takes some
    minutes or some hours for \~ 20,000 genes to obtain all genetic and
    gene-set parameters.

2.  Using genetic parameters from previous studies or from 1) to run
    `DECO`. It takes approximately 3 seconds for one gene-set. This
    approach can be used to analyze some thousands of gene-sets.

Run `DECO` for both variant and gene-set data.
----------------------------------------------

Below is an example for jointly analyzing variant data of schizophrenia
(SCZ) and a gene-set. The SCZ data is from
[https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-017-0497-y](Nguyen%20et%20al.,%202017,%20Genome%20Medicine).

### Data format

The data should be formatted as the `xInputGS` table below:

-   Mutation rates should start with `mut_`.

-   De novo mutations should start with `dn_`.

-   Case variant data should start with `cc_case`.

-   Control variant data should start with `cc_control`.

-   Gene-set data should be coded 1 (presence) or 0 (absence) and in the
    column `GS`.

``` r
xInputGS <- read.table("data/extTADA_SCZ_constrainedGenes.csv", header = TRUE, as.is = TRUE)
head(xInputGS)
```

### Use `DECO` function

The only function to obtain final results is `DECO`.

``` r
##SET PARAMETERS for DECO
nIteration = 1000 ###This should be higher (e.g, >=5,000)
nChain = 2
nCore = nChain
ntrio = c(1077, 1077, 1077) #Trio numbers: there are three categories of de novo mutations 
ncase = c(3157, 1091, 1353) #Case numbers: there are three population samples of rare case/control variants
ncontrol = c(4672, 1193, 4769) #Control numbers: there are three population samples of rare case/control variants

###RUN DECO
source('script/DECO.R') #Load the source code

###Obtain results
outSCZ <- DECO(
        inputData = xInputGS, ## Input data should be formatted as above
        Ndn = array(c(ntrio)), #rep(ntrio, 1), ##Two de novo categories
        Ncase = array(ncase), #rep(N$ca, 1), ##Two case categories
        Ncontrol = array(ncontrol), #rep(N$cn, 1), ##Two control categories
        nIteration = nIteration,# nIteration ## Number of iterations: should be upto higher for real data
        nChain = nChain, #Number of MCMC chains
        nCore = nCore #Number of computer cores
    )
```

### Outputs of `DECO`

Output is a list of:

1.  `dataPP`: all input data, posterior probabilities (PP) and FDRs of
    genes.

2.  `pars`: estimates of parameters.

Note: *π* is not estimated directly, but *α* was estimated with
*π* = *e*<sup>*α*</sup>/(*e*<sup>1 + *α*</sup>).

1.  `genesetInfo`: p-value (`pGS`), `q1` (q01) and `q0` (q00) of the
    tested gene-set.

2.  `mcmcData`: MCMC results of parameters.

We will take a look at these output information.

``` r
outSCZ$genesetInfo
fdr1 <- outSCZ$dataPP
dim(fdr1[fdr1$FDR < 0.05, ])
```

Using genetic parameters from previous studies or from 1) to run `DECO`
-----------------------------------------------------------------------

We describe an example in which we used genetic parameters from 1) to
run `DECO`. To use this approach, we need Bayes Factors (BF) and *π*.

``` r
outSCZ$pars
```

#### Read a gene-set (e.g. pLI09 in this example)

``` r
gsData <- readLines('data/pLI09.txt')
xInputGS2 <- fdr1[, 1:14]
gsInput <- rep(0, dim(xInputGS2)[1])
gsInput[is.element(xInputGS2[, 1], gsData)] <- 1

table(gsInput)

pi0 <- as.numeric(outSCZ$pars[1, 'Mode'])
xInputGS2[, 'GS'] <- gsInput
head(xInputGS2, 2)
```

#### Run `DECO`

``` r
outSCZ2 <- DECO(
  inputData = xInputGS2, ## Input data should be formatted as above
  pi0 = pi0,
  Ndn = array(c(ntrio)), #rep(ntrio, 1), ##Two de novo categories
  Ncase = array(ncase), #rep(N$ca, 1), ##Two case categories
  Ncontrol = array(ncontrol)
    )
```

#### Results

``` r
fdr2 <- outSCZ2$dataPP
fdr2[fdr2$FDR < 0.05, ]
outSCZ2$genesetInfo
```

``` r
stan_trace(outSCZ$mcmcData)
```

Citation
========

Tan-Hoang Nguyen, Xin He, Ruth Brown, Bradley T. Webb, Eli A. Stahl,
Kenneth S. Kendler, Vladimirov Vladimir, Brien P. Riley, Silviu-Alin
Bacanu. *DECO: Prioritizing risk genes for schizophrenia by jointly
analyzing rare variants and biological pathways*. Under Review (2020).
