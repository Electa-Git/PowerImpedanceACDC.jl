# System structure

```@meta
CurrentModule = PowerImpedanceACDC
```

## Network
The specified system is given in Network format and it consists of the components.
```@docs
PowerImpedanceACDC.Network
```

Network presents a struct with the collection of all components and their interconnections. Components are defined as type `Element` and the nodes present the dictionary, in which each node has its `Symbol` and it consists of the element pins connected to it.

`Element` is defined with the struct.
```@docs
PowerImpedanceACDC.Element
```

Network is initialized with the macro network as it follows:
```@docs
PowerImpedanceACDC.@network
```

### Add, delete elements and connect and disconnect element pins
```@docs
PowerImpedanceACDC.add!(n::Network, elem::Element)
```

```@docs
PowerImpedanceACDC.add!(n::Network, designator::Symbol, elem::Element)
```

```@docs
PowerImpedanceACDC.delete!(n::Network, designator::Symbol)
```

```@docs
PowerImpedanceACDC.connect!(n::Network, pins::Union{Symbol,Tuple{Symbol,Any}}...)
```

```@docs
PowerImpedanceACDC.disconnect!(n::Network, p::Tuple{Symbol,Symbol})
```



### Add, connect pins
```@docs
PowerImpedanceACDC.add!(n::Network, pin::Tuple{Symbol, Symbol})
```

```@docs
PowerImpedanceACDC.connect!(n::Network)
```

### Checking element connections
```@docs
PowerImpedanceACDC.check_lumped_elements(net :: Network)
```

## Form element from network
```@docs
PowerImpedanceACDC.composite_element
```

## Check port impedance and system stability

```@docs
PowerImpedanceACDC.determine_impedance
```
