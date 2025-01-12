[![linux-macos-windows](https://github.com/j-fu/LiquidElectrolytes.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/j-fu/LiquidElectrolytes.jl/actions/workflows/ci.yml)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://j-fu.github.io/LiquidElectrolytes.jl/dev)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://j-fu.github.io/LiquidElectrolytes.jl/stable)


LiquidElectrolytes.jl - A Generalized Poisson-Nernst-Planck solver
==================================================================

This package is a solver for generalized Poisson-Nernst-Planck models taking into account finite ion size, solvation and electroosmotic pressure based on formulations [derived from first principles of nonequilibrium thermodymamics](https://doi.org/10.1016/j.elecom.2014.03.015).  It utilizes a [thermodynamically consistent finite volume space discretization approach](https://doi.org/10.1007/s00211-022-01279-y) which uses the sum of the electrotstatic potential and the excess chemical potential as convective terms. It is realized on top of the solver kernel [VoronoiFVM.jl](https://github.com/WIAS-PDELib/VoronoiFVM.jl)  for coupled nonlinear PDE systems which takes advantage of [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) to generate the full linearization of the coupled nonlinear system as a basis for a robust Newton solver for the discrete nonlinear system.
