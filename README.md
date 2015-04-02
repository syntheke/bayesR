bayesR
======

Bayesian hierarchical model for complex trait analysis

### Quick start

#### Clone:

```sh
git clone https://github.com/syntheke/bayesR.git
```

#### Compile:

in the src folder
```sh
gfortran -o bayesR RandomDistributions.f90 baymods.f90 bayesR.f90
```

#### Run:
The example is from the [14th QTL-MAS workshop](http://jay.up.poznan.pl/qtlmas2010/index.html).
The simdata file in the example folder.
```sh
bayesR -bfile example/simdata -out simout -numit 10000 -burnin 5000 -seed 333
```

#### Help:

```sh
bayesR –help
```

#### Tell me more:

[BayesRmanual-0.75.pdf](https://github.com/syntheke/bayesR/blob/master/doc/BayesRmanual-0.75.pdf?raw=true)

#### New:
A more efficient implementation of the BayesR model for larger data sets is provided in the software bayesRv2. Source code to be released soon. For details see: [fast(er)BayesR.pdf](https://github.com/syntheke/bayesR/blob/master/doc/fast(er)BayesR.pdf?raw=true)

