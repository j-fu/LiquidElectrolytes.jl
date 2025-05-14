## General remarks
All physical quantities are assumed to be consistently represented through their values expressed in basic SI units
(m, kg, s, A, K, mol, cd), supported by the [LessUnitful.jl](https://j-fu.github.io/LessUnitful.jl/) package
built on top of [Unitful.jl](https://github.com/PainterQubits/Unitful.jl). 

## Electrolyte data
```@docs
AbstractElectrolyteData
ElectrolyteData
update_derived!
```

The default values for electrolyte data are those of an symmetric 0.1M aqueous binary electrolyte at 
298.5K with solvation number κ=10, ion molar volumes and masses similar to those of water molecules and
diffusion coefficients 2.0e-9 ``m^2/s``. All values given in SI base units:
```@example
using LiquidElectrolytes
ElectrolyteData()
```

```@docs
dlcap0(::ElectrolyteData)
debyelength(::ElectrolyteData)
chargedensity
chemical_potential
chemical_potentials!
rrate
iselectroneutral
isincompressible
c0_barc
``` 

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
### Upwind fluxes
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
