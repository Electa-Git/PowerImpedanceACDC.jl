# HVDCstability

Dev:
[![Build Status](https://travis-ci.com/Electa-Git/HVDCstability.jl.svg?token=wBsNbd12XnPoP4bx78Cy&branch=master)](https://travis-ci.com/Electa-Git/HVDCstability.jl)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://electa-git.github.io/HVDCstability.jl/latest/)

Tutorial about simulator construction and usage can be found on [link](https://github.com/Electa-Git/HVDCstability.jl/blob/master/HVDCstability.pdf).

## Download and install
To download the package, clone it using git. Afterwards just extract it and open the project inside julia.
The package can be added to all julia packages using the following command.
```
using Pkg
Pkg.add(PackageSpec(path = pwd()))
```
In order to run the package it is necessary to install supporting packages: `Parameters`, `DataStructures`
`FileIO`, `DelimitedFiles`, `Plots`, `LaTeXStrings`, `DSP`, `Compat`, `SymEngine`, `LinearAlgebra`, `NLsolve`, `ForwardDiff`. Each of them can be added by calling:
```
using Pkg
Pkg.add("package_name")
```

To remove the package, call the following command:
```
using Pkg
Pkg.rm("HVDCstability")