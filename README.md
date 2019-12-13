# HVDCstability

Dev:
[![Build Status](https://api.travis-ci.com/Aleksandra-Lekic/HVDCstability.jl.svg?token=8MGGs8A1RyNuFsSbtpur&branch=master)](https://travis-ci.com/Aleksandra-Lekic/HVDCstability.jl)

Tutorial about simulator construction and usage can be found on [link](https://github.com/Aleksandra-Lekic/HVDCstability.jl/blob/master/HVDCstability.pdf).

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
```
