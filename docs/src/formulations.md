# System structure

```@meta
CurrentModule = HVDCstability
```

## Network
The specified system is given in Network format and it consists of the components.
```@docs
HVDCstability.Network
```

Network presents a struct with the collection of all components and their interconnections. Components are defined as type `Element` and the nodes present the dictionary, in which each node has its `Symbol` and it consists of the element pins connected to it.

`Element` is defined with the struct.
```@docs
HVDCstability.Element
```

Network is initialized with the macro network as it follows:
```@docs
HVDCstability.@network
```

### Add, delete elements and connect and disconnect element pins
```@docs
HVDCstability.add!(n::Network, elem::Element)
```

```@docs
HVDCstability.add!(n::Network, designator::Symbol, elem::Element)
```

```@docs
HVDCstability.delete!(n::Network, designator::Symbol)
```

```@docs
HVDCstability.connect!(n::Network, pins::Union{Symbol,Tuple{Symbol,Any}}...)
```

```@docs
HVDCstability.disconnect!(n::Network, p::Tuple{Symbol,Symbol})
```



### Add, connect pins
```@docs
HVDCstability.add!(n::Network, pin::Tuple{Symbol, Symbol})
```

```@docs
HVDCstability.connect!(n::Network)
```

### Checking element connections
```@docs
HVDCstability.check_lumped_elements(net :: Network)
```

## Form element from network
```@docs
HVDCstability.composite_element
```

## Check port impedance and system stability

```@docs
HVDCstability.determine_impedance
```
