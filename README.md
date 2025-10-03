## Installation

The package can be installed directly from github, using either `devtools` or `remotes`:

```R
devtools::install_github("https://github.com/HervePerdry/AMVCTpaper")
```
or 
```R
remotes::install_github("https://github.com/HervePerdry/AMVCTpaper")
```

It can then be loaded with 

``` r
library(AMVCTpaper)
```

### Overview

This package contains a set of functions to calculate the evolution and equilibrium values of quantities of interest
under the AMVTC model (co-occurrence of Assortative Mating and Vertical Cultural Transmission) as described here:
https://www.biorxiv.org/content/10.1101/2023.04.08.536101v3. This includes both theoretical calculations as well as the
simple forward-time simulations detailed in the publication; as well as code to reproduce the main figures therein.

#### Theoretical evolution equations

The evolution of all parameters in the AMVTC model can be estimated with the function pop.evolution() for  a population
which is initially at gametic equilibrium, without assortative mating. The function takes the following arguments

|          | |
|----------|--------------------------------------------------------------------------------------------|
| `g0`     | Standard deviation of gametic value in a population without assortative mating  |
| `e`      | Standard deviation of environmental effects  |
| `r.ho`   | Correlation between mates   |
| `nu`     | Correlation between parents' and offspring environmental effects  |
| `N`      | Number of causal variants, all will have a minor-allele frequency of 0.5 and are unlinked  |
| `nb.gen` | Number of generations  |


The function returns a data frame with class "pop.evolution", with columns 

|           | |
|-----------| ---------------------------------------------------|
| `t`       | the generation, from 0 to nb.gen                   |
| `N.kappa` | the value of $N\overline\kappa(t)$ |
| `g`       | the standard deviation of the gametic value $g(t)$ |
| `r.ga`    | the gametic correlation, $r_ga(t)$   |
| `rho`     | the gene-environment correlation $\rho(t)$ |
| `a`       | the standard deviation of the genetic value $a(t)$ |
| `sigma2`  | the variance of the phenotype, $\sigma^2(t)$ |

Example:


``` r
pop.evolution(g0 = sqrt(0.5), e = 1, r.ho = 0.6, nu = 0.4, N = 1e4, nb.gen = 10)
```

```
## Evolution of a population with parameters :
##           g0            e         r.ho           nu            N 
## 7.071068e-01 1.000000e+00 6.000000e-01 4.000000e-01 1.000000e+04 
## 
##     t  N.kappa         g      r.ga        rho        a   sigma2
## 1   0 1.000000 0.7071068 0.0000000 0.00000000 1.000000 2.000000
## 2   1 1.000000 0.7071068 0.1500000 0.08607737 1.072381 2.334615
## 3   2 1.074993 0.7331414 0.1844841 0.11965463 1.128411 2.543351
## 4   3 1.136646 0.7538719 0.2058256 0.13419146 1.170725 2.684799
## 5   4 1.185287 0.7698333 0.2200194 0.14136424 1.202528 2.786061
## 6   5 1.223023 0.7819921 0.2299400 0.14538415 1.226477 2.860866
## 7   6 1.252109 0.7912359 0.2370940 0.14788161 1.244579 2.917077
## 8   7 1.274473 0.7982710 0.2423556 0.14955031 1.258312 2.959711
## 9   8 1.291659 0.8036352 0.2462754 0.15071946 1.268764 2.992218
## 10  9 1.304866 0.8077331 0.2492212 0.15156363 1.276740 3.017081
## 11 10 1.315017 0.8108689 0.2514488 0.15218498 1.282839 3.036134
```

#### Equilibrium values

While `pop.evolution()` gives the theoretical evolution of all parameters but not the equilibrium values. These are
provided by the function `AMVCT()`, which calculated the theoretical equilibrium values of all parameters in the model, as
well as the resulting correlations between various pairs of genetic and non-genetic components for members of a trio;
described in Tables 1 and 2 in the publication. The function will solve the non-linear system of equations with the
function `solve.a.rho()`. This function has arguments `h2.0`, `g0`, `e`, `r.ho`, and `nu`. The
argument `h2.0` allow to give the value of the heritability in the initial population, without AM and VCT.
All other arguments are similar to the arguments of `pop.evolution`.
The user must provide either `g0` and `e`, or `h2.0`, but not both.

`AMVCT()` returns an object of class `AMVCT`, which is a list with components g0`, `e`, `r.ho`, `nu`, `rho`, `a, `sigma2`,
`r.ga`, `g`, `h2.SNP`, `decompose.sigma2`, `cor.mates` (corresponding to table 1 of the supplementary material)
`cor.parent.offspring` (corresponding to table 2 of the supplementary material) `cor.pheno.parent.offspring` (correlation
of the phenotypes of parent and offspring).

Example:

``` r
AMVCT(g0 = sqrt(0.5), e = 1, r.ho = 0.6, nu = 0.4)
```

```
## Assortative mating and verticulal cultural transmission with parameters
## g0 = 0.7071068 e =  1 r.ho =  0.6 nu = 0.4 
## Corresponding to an heritability without AM and VCT h²_0 = 0.5 
## ----------------------------------------------------------------
## Values at equilibrium:
## Gametic correlation r_ga =  0.2586618 
## Gametic effect g =  0.8212527 or g^2 = 0.674456 
## Total additive effect a = 1.302989 or a^2 =  1.69778 
## Correlation (A, E) rho = 0.1541091 
## Phenotype variance sigma^2 = 3.099385 decomposing as
##             a2.a deux.rho.a.e.rho               e2 
##        1.6977796        0.4016049        1.0000000 
## 
## SNP-heritability h^2_SNP = 0.685018 
## 
## Correlations between mates
##           A1        E1        A2        E2
## A1 1.0000000 0.1541091 0.4110108 0.3387163
## E1 0.1541091 1.0000000 0.3387163 0.2791380
## A2 0.4110108 0.3387163 1.0000000 0.1541091
## E2 0.3387163 0.2791380 0.1541091 1.0000000
## 
## Correlations between parent and offspring
##           A3        E3
## A1 0.7055054 0.1541091
## E1 0.2464127 0.4000000
## 
## Leading to a correlation between parent's and offspring phenotypes = 0.6838997
```

#### Forward time simulations

The evolution of all parameters can also be estimated from a simple simulated population using pop.sim(), with a fixed
size across generations, with both AM and VCT starting in generation zero. Here, we apply a simple model for forming M/2
couples from the population of M individuals, each of which produce exactly 2 offspring in order to form the next
generation; with all individuals having a simulated phenotype value following the AMVCT model. The function takes the
the same arguments as `pop.evolution`, and a few new arguments:

|              | |
|--------------|---------------------|
| pop.size     | Population size  |
| digest       | Logical. Set to `TRUE` to send back only a digest of the results  |
| keep.N.kappa | Logical. Set to `TRUE` to keep trace of $N \overline\kappa(t)$.  |

Beware: setting `keep.N.kappa` to `TRUE` results in lengthy computations.

This function will return, if `digest` is `TRUE`, a data frame similar to the result of `pop.evolution`, with an
additional column for `e` which contain the sd of the environmental components at each generation (this fluctuates
slightly around the value of the corresponding parameter, due to random sampling). If `digest` is `FALSE`, it returns a list with
member `digest` (the digest as described before), `G` (the matrix of genotypes at the last generation), `A` (the genetic
value at the last generation), `E` (the environmental value at the last generation), `P` (the phenotypes at the last
generation).

Example of a forward-time simulation over 20 generations, with N = 100 causal loci and a population of size 25,000:


``` r
R <- pop.sim(g0 = sqrt(0.5), e = 1, r.ho = 0.6, nu = 0.4, pop.size = 25000, N = 100, nb.gen = 20, TRUE, TRUE)
R
```

```
## Evolution of a population with parameters :
##           g0            e         r.ho           nu            N     pop.size 
## 7.071068e-01 1.000000e+00 6.000000e-01 4.000000e-01 1.000000e+02 2.500000e+04 
## 
##     t         e  N.kappa         g        r.ga         rho         a   sigma2
## 1   0 1.0002658 1.004605 0.7069211 -0.00511479 0.001421885 0.9997374 2.002850
## 2   1 1.0045939 1.012339 0.7104319  0.14767821 0.089338428 1.0762229 2.360644
## 3   2 1.0005464 1.082731 0.7359697  0.17833003 0.123079808 1.1292725 2.554483
## 4   3 0.9999503 1.140073 0.7573535  0.22156733 0.139767915 1.1752395 2.709594
## 5   4 0.9988963 1.188276 0.7678751  0.21473912 0.144987184 1.2054171 2.799979
## 6   5 1.0003705 1.229643 0.7852659  0.22053569 0.140389921 1.2271906 2.851436
## 7   6 0.9981794 1.247136 0.7903801  0.23305861 0.146285678 1.2414236 2.900038
## 8   7 0.9995990 1.274326 0.7984058  0.23954164 0.146502473 1.2564527 2.945871
## 9   8 0.9889585 1.279561 0.7991638  0.24836570 0.142602857 1.2602772 2.921807
## 10  9 0.9993944 1.296748 0.8041249  0.24788938 0.150539179 1.2732699 3.003127
## 11 10 1.0006085 1.305455 0.8098158  0.26831286 0.151686817 1.2841373 3.040036
## 12 11 0.9998658 1.324382 0.8137290  0.24771172 0.149554744 1.2867232 3.040208
## 13 12 1.0054871 1.328762 0.8121542  0.24093768 0.156526541 1.2862052 3.070188
## 14 13 1.0022211 1.325847 0.8144139  0.26541649 0.157304739 1.2965505 3.094303
## 15 14 1.0032238 1.343101 0.8202577  0.25241766 0.153079123 1.2966834 3.086116
## 16 15 1.0014692 1.340683 0.8172408  0.26525332 0.149530449 1.3021828 3.088624
## 17 16 1.0002602 1.344086 0.8192536  0.26432510 0.156922359 1.3025240 3.105986
## 18 17 0.9933717 1.344175 0.8180518  0.26448349 0.158460989 1.3039709 3.097645
## 19 18 1.0028339 1.362602 0.8233094  0.26579733 0.161343035 1.3096247 3.144588
## 20 19 1.0001487 1.358410 0.8228374  0.25824976 0.152866773 1.3070760 3.108422
## 21 20 1.0012519 1.349128 0.8179382  0.25924757 0.156191479 1.3025573 3.106567
```

We can plot a comparison to the theoretical evolution values of `pop.evolution()` and equilibrium values of `AMVCT()`:

``` r
ev <- pop.evolution(g0 = sqrt(0.5), e = 1, r.ho = 0.6, nu = 0.4, N = 100, nb.gen = 20) 
limits <- AMVCT(g0 = sqrt(0.5), e = 1, r.ho = 0.6, nu = 0.4)

par(mfrow=c(2,2), cex = .7, mai = c(3,3,1,1)/5)
## plotting evolution of rho
plot(R$t, R$rho, type = "l", xlab = "t", ylab = expression(rho))
lines(ev$t, ev$rho, col = "red")
abline(h = limits$rho, col = "red", lty = 3)

## plotting evolution of N kappa bar
plot(R$t, R$N.kappa, type = "o", xlab = "t", ylab = expression(N * bar(kappa)))
lines(ev$t, ev$N.kappa, col = "red")
abline(h = 1/(1-limits$r.ga), col = "red", lty = 3)

## plotting evolution of a
plot(R$t, R$a, type = "o", xlab = "t", ylab = "a") 
lines(ev$t, ev$a, col = "red")
abline(h = limits$a, col = "red", lty = 3)

## plotting evolution of r.ga
plot(R$t, R$r.ga, type = "o", xlab = "t", ylab = expression(r[ga]))
lines(ev$t, ev$r.ga, col =" red")
abline(h = limits$r.ga, col = "red", lty = 3)
```

![plot](./figure/plots-1.png)

#### Recreating the plots in the article

The directory `figures/` contains instruction and script fpr
reproducing all figures in the paper.

Running `simus-50.r` will re-generate the files `simus-50-N1000-pop100k.rds` and `simus-50-N1000-pop25k.rds` 
that are needed for some of the figures.

Running `figures-AMVCT-paper.r` will generate the figures.

The files `multiple.pop.sim.r`, `fig.evol.r` and `fig.nu.r` contain functions definition that 
are needed by `simus-50.r` and `figures-AMVCT-paper.r`.
