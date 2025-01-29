# Impedance & Stability

PowerImpedanceACDC enables saving the node impedance and components data as plots and
in the text file.

```@meta
CurrentModule = PowerImpedanceACDC
```
## Determine impedance

## Saving data
Component data can be saved in text files using
```@docs
PowerImpedanceACDC.save_data
```



## Plotting options

### Component data
Component specific data can be plotted using
```@docs
PowerImpedanceACDC.plot_data
```

### Bode plot
I think second option is the newer version with added functionalities of Julia plotting package, so first one can be deleted.
```@docs
PowerImpedanceACDC.bode
```

```@docs
PowerImpedanceACDC.bodeplot
```

### Nyquist plot

```@docs
PowerImpedanceACDC.nyquistplot
```

```@docs
PowerImpedanceACDC.EVD
```

```@docs
PowerImpedanceACDC.passivity
```

```@docs
PowerImpedanceACDC.small_gain
```