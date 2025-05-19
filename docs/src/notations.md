# Notations

## Remark on Units
All physical quantities are assumed to be consistently represented through their values expressed in basic SI units
(m, kg, s, A, K, mol, cd), supported by the [LessUnitful.jl](https://j-fu.github.io/LessUnitful.jl/) package
built on top of [Unitful.jl](https://github.com/PainterQubits/Unitful.jl). Values of physical constants are 
obtained via [LessUnitful.jl](https://j-fu.github.io/LessUnitful.jl/) from 
[PhyscialConstants.CODATA2018](https://github.com/JuliaPhysics/PhysicalConstants.jl).


## Basic unknowns
`u` refers to a solution vector,
`electrolyte` refers to an instance of [`ElectrolyteData`](@ref)

| Symbol                                 | Unit                        | Meaning                 | Notation in code                          |
|:---------------------------------------|:----------------------------|:------------------------|:------------------------------------------|
| ``\phi ``                              | ``V``                       | electrostatic potential | `u[iϕ]`                                   |
| ``p``                                  | ``Pa``                      | pressure                | `u[ip]`                                   |
| ``c_i``                                | ``\frac{mol}{m^3}``         | molar concentration     | `u[i]`                                    |
| ``\gamma_i=\gamma_i(c_1\dots c_n, p)`` | ``1``                       | activity coefficient    | see `electrolyte.actcoeff!`               |
| ``N_i``                                | ``\frac{mol}{m^2\cdot{}s}`` | molar flux              | see [`LiquidElectrolytes.pnpflux!`](@ref) |
| ``q=F\sum_{i=1}^N z_i c_i``            | ``\frac{A\cdot{}s}{m^3}``   | electric charge density | see [`chargedensity`](@ref)               |

## Species & Mixture data
`electrolyte` refers to an instance of [`ElectrolyteData`](@ref)


| Symbol                    | Unit                | Meaning                         | Notation in code                     |
|:--------------------------|:--------------------|:--------------------------------|:-------------------------------------|
| ``z_i``                   | 1                   | charge number                   | `electrolyte.z[i]`                   |
| ``D_i``                   | ``\frac{m^2}{s}``   | diffusion coefficient           | `electrolyte.D[i]`                   |
| ``v_i``                   | ``\frac{m^3}{mol}`` | molar volume                    | `electrolyte.v0`, `electrolyte.v[i]` |
| ``M_i``                   | ``\frac{kg}{mol}``  | molar mass                      | `electrolyte.M0`, `electrolyte.M[i]` |
| ``κ_i``                   | ``1``               | solvation number                | `electrolyte.κ[i]`                   |
| ``ε``                     | ``1``               | mixture dielectric permittivity | `electrolyte.ε`                      |
| ``T``                     | ``K``               | temperature                     | `electrolyte.T`                      |
| ``c^\circ=\frac{1}{v_0}`` | ``\frac{mol}{m^3}`` | reference concentration         | `1/electrolyte.v0`                   |
| ``p^\circ``               | ``Pa``              | reference pressure              | `electrolyte.p_bulk`                 |
| ``\mu_i^\circ``           | ``\frac{J}{mol}``   | reference chemical potential    | 0                                    |


## Physical constants
`electrolyte` refers to an instance of [`ElectrolyteData`](@ref)

| Symbol                         | Unit                            | Meaning                        | Notation in code |
|:-------------------------------|:--------------------------------|:-------------------------------|:-----------------|
| ``ε_0=8.854187\cdot 10^{-12}`` | ``\frac{A\cdot{}s}{V\cdot{}m}`` | vacuum dielectric permittivity | `electrolyte.ε0` |
| ``F=9.64853321 \cdot 10^4``    | ``\frac{A\cdot{}s}{mol}``       | Faraday constant               | `electrolyte.F`  |
| ``R=8.31446261815324``         | ``\frac{J}{mol\cdot{}K}``       | molar gas constant             |                  |

## Derived data
`electrolyte` refers to an instance of [`ElectrolyteData`](@ref)


| Symbol                                                                                                | Unit                | Meaning                        | Notation in code      |
|:------------------------------------------------------------------------------------------------------|:--------------------|:-------------------------------|:----------------------|
| ``\bar v_i = v_i + \kappa_i v_0``                                                                     | ``\frac{m^3}{mol}`` | solvated ion molar volume      | `electrolyte.barv[i]` |
| ``c_0 = \frac{1}{v_0}- \sum_{i=1}^N \frac{\bar v_i}{v_0}c_i``                                         | ``\frac{mol}{m^3}`` | molar concentration of solvent | see [`c0_barc`](@ref)  |
| ``\bar c = \sum_{i=0}^N  c_i= \frac{1}{v_0}+ \sum_{i=1}^N \left (1 -\frac{\bar v_i}{v_0}\right) c_i`` | ``\frac{mol}{m^3}`` | molar concentration of mixture | see [`c0_barc`](@ref)  |
| ``m_i=\frac{M_i+\kappa_iM_0}{M_0}=\frac{M_i}{M_0}+\kappa_i``                                          | ``\frac{kg}{mol}``  | solvated molar mass ratio      | `electrolyte.Mrel`    |


    
