# HVDCstability.jl Documentation

```@meta
CurrentModule = HVDCstability
```

## Overview

HVDCstability.jl is a Julia package, which focuses on the impedance-based stability. Simulator provides support for the power system construction using detailed components. Components and their interconnections are modeled as multiport networks using ABCD parameters.

Implemented components are:
- Impedances
- AC and DC grid equivalents (sources)
- Transformers
- Overhead lines
- Cables
- MMCs (using various controllers)

Developed by:
- Aleksandra Lekic, Philippe De Rua, Jef Beerten - KU Leuven / EnergyVille
- with the help of:
                    - Hakan Ergun and Jay Dave for the power flow implementation

                    - Thomas Roose, Ozgur Can Sakinci for the MMC and MMC controls modeling

                    - Willem Leterme for the overhead line and autotransformer model

                    - Alejandro Bajo Salas, whose thesis presents the initial idea for this simulator

## Special Thanks To
This work is part of the Neptune project, supported by the Energy Transition Fund, FOD Economy, Belgium.  
