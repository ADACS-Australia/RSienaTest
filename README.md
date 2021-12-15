# MPI-optimised version of RSienaTest

Original source from https://rdrr.io/rforge/RSienaTest/, based on version 1.2-30 (July 13, 2021).

There is an install_github function to install R packages hosted on GitHub in the `devtools` package.

```
install_github("ADACS-Australia/RSienaTest")
```

Needs [mpi-Rscript](https://github.com/ADACS-Australia/mpi-Rscript) to run properly.

Example usage:
```
mpirun -n 4 mpi-Rscript example.R
```

example.R:
```
library(RSienaTest)

# load data and initial estimation
load('example.RData')

groupModel.ec <- sienaBayes(GroupsModel, data = example,
                            effects = GroupEffects, priorMu = Mu, priorSigma = Sig,
                            priorKappa = 0.01,
                            prevBayes = groupModel.e,
                            nmain=3, nrunMHBatches=40,
                            silentstart=FALSE,
                            clusterType="MPI")
```
