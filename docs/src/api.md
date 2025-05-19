# Model implementations

## Poisson-Nernst-Planck system

```@docs
AbstractElectrochemicalSystem
PNPSystem
LiquidElectrolytes.DGL_gamma!
electrolytedata
solventconcentration
chemical_potentials
electrochemical_potentials
```
## Upwind fluxes
```@docs
LiquidElectrolytes.μex_flux!
LiquidElectrolytes.act_flux!
LiquidElectrolytes.cent_flux!
```

## Poisson-Boltzmann system
```@docs
PBSystem
```

## Poisson-Nernst-Planck-Stokes system
```@docs
PNPStokesSolver
```

## Utilities
```@docs
RLog
RExp
```
