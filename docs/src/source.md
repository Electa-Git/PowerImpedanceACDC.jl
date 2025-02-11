# Source models
```@meta
CurrentModule = PowerImpedanceACDC
```
Source models can be used to represent constant voltage loads or to construct Thevenin-equivalents.


## AC source
```@docs
PowerImpedanceACDC.ac_source
```
## DC source
A DC source has currently not been implemented. If a constant DC voltage is required, connect DC side of converter to MMC with DC voltage control.