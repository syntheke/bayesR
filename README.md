bayesR
======

Bayesian hierarchical model for complex trait analysis

(https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004969)

# Update 01/04/2021
Various new features and improvements:

* further reduced memory requirements 
* inclusion of covariates
* grouped effects models to fit more complex models (e.g. partioning of variance)
* flat input files (not supported in bayesRv2)
* fitted values
* prediction of phenotypes
* unified source code


Previous release moved to folder /old

### Quick start

#### Clone:

```sh
git clone https://github.com/syntheke/bayesR.git
```

#### Compile:

in the src folder
```sh
gfortran -o bayesR -O2 -cpp RandomDistributions.f90 baymods.f90 bayesR.f90
gfortran -o bayesRv2 -O2 -cpp –Dblock -fopenmp RandomDistributions.f90 baymods.f90 bayesR.f90
ifort -o bayesR -O3 -fpp RandomDistributions.f90 baymods.f90 bayesR.f90
ifort -o bayesRv2 -O3 -fpp -Dblock -openmp -static RandomDistributions.f90 baymods.f90 bayesR.f90
```

#### Run:

```sh
bayesR -bfile simdata -out simout
```
##### Example1
Example from the [14th QTL-MAS workshop](http://jay.up.poznan.pl/qtlmas2010/index.html).
```sh
bayesR -bfile example/simdata -out simout -numit 10000 -burnin 5000 -seed 333
```
##### Example2
Genome position specific priors
```sh
bayesR -bfile simdata2 -out simout2 -numit 10000 -burnin 5000 -seed 333 -n 2 -snpmodel mod2 -segment seg
```
##### Example3
Grouped effects with mixture priors
```sh
bayesR -bfile simdata2 -out simout3 -numit 10000 -burnin 5000 -seed 333 -n 2 -snpmodel mod3 -segments seg -varcomp var3
```

#### Help:

```sh
bayesR -help
```

#### Tell me more:

[BayesRManual.pdf](https://github.com/syntheke/bayesR/blob/master/doc/BayesRManual.pdf?raw=true)

