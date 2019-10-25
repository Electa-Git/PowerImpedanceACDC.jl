# HVDCstability.jl Documentation

```@meta
CurrentModule = HVDCstability
```

## Overview

HVDCstability.jl is a Julia package, which focuses on the impedance-based stability. Simulator provides support for the power system construction using detailed components. Components and their interconnections are modeled as multiport networks using ABCD parameters.

Implemented components are:
- Impedance
- AC and DC grid equivalents (sources)
- Transfromers
- Overhead lines
- Cables
- MMC

Developed by:
- Aleksandra Lekic, Jef Beerten - KU Leuven / EnergyVille
- with the help of: Willem Leterme, Alejandro Bajo Salas, Thomas Roose and Ozgur Can Sakici - KU Leuven / EnergyVille


## Installation of HVDCstability

The latest stable release of HVDCstability can be installed using the Julia package manager with

```julia
Pkg.clone("")
```

## Special Thanks To
This work is part of the Neptune project, supported by the Energy Transition Fund, FOD Economy, Belgium.  
