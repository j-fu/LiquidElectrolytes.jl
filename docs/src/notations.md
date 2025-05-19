# Notations

## Remark on Units
All physical quantities are assumed to be consistently represented through their values expressed in basic SI units
(m, kg, s, A, K, mol, cd), supported by the [LessUnitful.jl](https://j-fu.github.io/LessUnitful.jl/) package
built on top of [Unitful.jl](https://github.com/PainterQubits/Unitful.jl). Values of physical constants are 
obtained via [LessUnitful.jl](https://j-fu.github.io/LessUnitful.jl/) from 
[PhyscialConstants.CODATA2018](https://github.com/JuliaPhysics/PhysicalConstants.jl).


## Basic unknowns

| Symbol                                                       | Unit                            | Meaning                         |
|:-------------------------------------------------------------|:--------------------------------|:--------------------------------|
| ``\phi ``                                                    | ``V``                           | electrostatic potential         |
| ``p``                                                        | ``Pa``                          | pressure                        |
| ``c_i``                                                      | ``\frac{mol}{m^3}``             | molar concentration             |
| ``\gamma_i=\gamma_i(c_1\dots c_n, p)``                       | ``1``                           | activity coefficient            |
| ``N_i``                                                      | ``\frac{mol}{m^2\cdot{}s}``     | molar flux                      |
| ``q=F\sum_{i=1}^N z_i c_i``                                  | ``\frac{A\cdot{}s}{m^3}``       | electric charge density         |

## Species & Mixture data
| Symbol                                                       | Unit                            | Meaning                         |
|:-------------------------------------------------------------|:--------------------------------|:--------------------------------|
| ``z_i``                                                      | 1                               | charge number                   |
| ``D_i``                                                      | ``\frac{m^2}{s}``               | diffusion coefficient           |
| ``v_i``                                                      | ``\frac{m^3}{mol}``             | molar volume                    |
| ``M_i``                                                      | ``\frac{kg}{mol}``              | molar mass                      |
| ``\kappa_i``                                                 | ``1``                           | solvation number                |
| ``\varepsilon``                                              | ``1``                           | mixture dielectric permittivity |
| ``T``                                                        | ``K``                           | temperature                     |
| ``c^\circ=\frac{1}{v_0}``                                    | ``\frac{mol}{m^3}``             | reference concentration         |
| ``p^\circ``                                                  | ``Pa``                          | reference pressure              |
| ``\mu_i^\circ``                                              | ``\frac{J}{mol}``               | reference chemical potential    |

## Physical constants
| Symbol                                                       | Unit                            | Meaning                         |
|:-------------------------------------------------------------|:--------------------------------|:--------------------------------|
| ``\varepsilon_0=8.8541878188\cdot 10^{-12}``                 | ``\frac{A\cdot{}s}{V\cdot{}m}`` | vacuum dielectric permittivity  |
| ``F=9.64853321 \cdot 10^4``                                  | ``\frac{A\cdot{}s}{mol}``       | Faraday constant                |
| ``R=8.31446261815324``                                       | ``\frac{J}{mol\cdot{}K}``       | molar gas constant              |

## Derived data
| Symbol                                                                                                | Unit                | Meaning                        |
|:------------------------------------------------------------------------------------------------------|:--------------------|:-------------------------------|
| ``\bar v_i = v_i + \kappa_i v_0``                                                                     | ``\frac{m^3}{mol}`` | solvated ion molar volume      |
| ``c_0 = \frac{1}{v_0}- \sum_{i=1}^N \frac{\bar v_i}{v_0}c_i``                                         | ``\frac{mol}{m^3}`` | molar concentration of solvent |
| ``\bar c = \sum_{i=0}^N  c_i= \frac{1}{v_0}+ \sum_{i=1}^N \left (1 -\frac{\bar v_i}{v_0}\right) c_i`` | ``\frac{mol}{m^3}`` | molar concentration of mixture |
| ``m_i=\frac{M_i+\kappa_iM_0}{M_0}=\frac{M_i}{M_0}+\kappa_i``                                          | ``\frac{kg}{mol}``  | solvated molar mass ratio      |


    
