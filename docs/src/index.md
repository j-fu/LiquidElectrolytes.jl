# LiquidElectrolytes.jl

```@docs
LiquidElectrolytes
```

## Poisson-Nernst-Planck (PNP) system
The package is devoted to the numerical solution of the generalized Poisson-Nernst-Planck (PNP) system in the
following formulation:

Poisson equation:
```math
-\nabla \varepsilon\varepsilon_0 \nabla \phi = q = F\sum_{i=1}^N z_i c_i
```
Momentum balance taking derived from the Navier-Stokes equation in mechanical equilibrium:
```math
    -\Delta p + q\nabla \phi = 0
```
Continuity equation:
```math
\partial_t c_i  + \nabla \cdot \mathbf N_i  = 0  \qquad (i=1\dots N)
```
Generalized Nernst-Planck flux:
```math
 N_i=-D_i c_i \left( \nabla\log \left(\gamma_i\frac{c_i}{c^\circ}\right) + z_i\frac{F}{RT}\nabla\phi\right)  \qquad  (i=1\dots N)
```
The activity coefficients ``\gamma_i = \gamma_i(c_1\dots c_N, p)`` encode a particular electrolyte model. In general,
    the fluxes ``N_i`` are coupled due to the dependency of ``\gamma_i`` on all concentrations and the pressure.
E.g. ``\gamma_i=1`` describes the classical Guoy-Chapman model. The default implements the model of
Dreyer, Guhlke, Müller, Landstorfer. See the section on [Activity coefficients](@ref) 
for more.

## Poisson-Boltzmann (PB) system
In thermodynamic equilibrium, the fluxes ``N_i`` of the PNP system are zero, which allows to replace the continuity
and flux equations by the algebraic condition
```math
c_i=c_i^b \frac{\gamma_i(c_1^b\dots c_N^b, p^b)}{\gamma_i(c_1\dots c_N,p)}\exp\left(z_i\frac{F}{RT}(\phi^b-\phi)\right),

```
where ``\phi^b, p^b, c_1^b\dots c_N^b`` are reference values usually imposed as boundary conditions at the interface
to the bulk of the electrolyte. Poisson equation and momentum balance are kept.
In general, this algebraic condition for given ``\phi, p`` comprises a nonlinear system of
equations due to the dependency of ``\gamma_i`` on all concentrations and the pressure.



