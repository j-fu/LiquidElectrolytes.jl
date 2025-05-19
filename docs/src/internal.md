# Internal API

## Electrolyte data


```@docs
LiquidElectrolytes.solvepressure
LiquidElectrolytes.edgevelocity
LiquidElectrolytes.value_or_getindex
LiquidElectrolytes.chargedensity!
LiquidElectrolytes.wnorm
```


## Poisson-Nernst-Planck
```@docs
LiquidElectrolytes.pnpreaction!
LiquidElectrolytes.pnpflux!
LiquidElectrolytes.pnpstorage!
LiquidElectrolytes.pnpbstorage!
LiquidElectrolytes.dμex
```
## Poisson-Boltzmann
```@docs
LiquidElectrolytes.pbreaction!
LiquidElectrolytes.pbflux!
```

## Poisson-Nernst-Planck-Stokes
```
LiquidElectrolytes.flowsolver
```

## Electrochemical calculations
```@docs
LiquidElectrolytes.splitz
```

## General tools
```@docs
LiquidElectrolytes.showstruct
LiquidElectrolytes.myround
```
