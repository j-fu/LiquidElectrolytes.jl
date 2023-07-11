# Changes
## Pre-Registration (versions <0.1.0)

### v0.0.23  BREAKING July 11, 2023
- Breaking: 
  - Fixed PNP discretization: for consistency with equilibrium we need to use  ``\bar c \nabla \frac{c_i}{\bar c} \dots``
    instead of ``\nabla c \dots``. 
  - [`ivsweep`](@ref)  and [`dlcapsweep`](@ref)  now have their own return types. These
    can be amended by more information without breaking codes
  - [`ivsweep`](@ref)  and [`dlcapsweep`](@ref)  have now a `store_solutions` keyword argument (false by default)
    which indicates if all solutions should be stored.
- Introduced Poisson-Boltzmann solver
- Surface species example (thanks @smaasz!)
- Checked consistency with equilibrium for all variants (equilibrium, pb, pnp)


### v0.0.22  BREAKING June 28, 2023
- Update return values in [`ivsweep`](@ref) 

### v0.0.21  BREAKING June 26, 2023
- More general return values in [`ivsweep`](@ref) 

