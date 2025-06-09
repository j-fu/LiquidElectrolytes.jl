"""
$(TYPEDEF)

Abstract super type for electrochemical systems
"""
abstract type AbstractElectrochemicalSystem end

"""
$(TYPEDEF)

Abstract super type for electrolytes.
"""
abstract type AbstractElectrolyteData end

"""
    DGML_gamma!(γ, c, p, electrolyte)

Activity coefficients according to Dreyer, Guhlke, Müller, Landstorfer.

```math
γ_i = \\exp\\left(\\frac{\\tilde v_i p}{RT}\\right) \\left(\\frac{\\bar c}{c_0}\\right)^{m_i} \\frac{1}{v_0\\bar c}
```

Input:
- c: vector of concentrations
- p: pressure
- electrolyte: instance of [`ElectrolyteData`](@ref)

Output: 
- γ (mutated) activity coefficients
"""
function DGML_gamma!(γ, c, p, electrolyte)
    (; Mrel, tildev, v0, RT, v0, cspecies, rexp) = electrolyte
    c0, barc = c0_barc(c, electrolyte)
    for ic in cspecies
        γ[ic] = rexp(tildev[ic] * p / RT) * (barc / c0)^Mrel[ic] * (1 / (v0 * barc))
    end
    return nothing
end


"""
$(TYPEDEF)

Data for electrolyte. It is defined using [`Base.@kwdef`](https://docs.julialang.org/en/v1/base/base/#Base.@kwdef)
allowing for keyword constructors like
```julia
    ElectrolyteData(nc=3,z=[-1,2,1])
```
The struct has three groups of fields:
- Problem parameters: these are meant to be set by user in order to provide simulation parameters.
- Derived values: these are derived from problem parameters by default, or via [`update_derived!`](@ref) after setting some of the problems pararmeters
- Fixed values: these are mostly physical constants and values derived from them, they should not be changed.
- Reserved fields: they are used to pass information during simulations and should not be set by the user.

Fields (reserved fields are modified by some algorithms):
$(TYPEDFIELDS)
"""
@kwdef mutable struct ElectrolyteData{Tγ, Tcache, Texp, Tlog, Tflux} <: AbstractElectrolyteData
    "Number of charged species ``N``."
    nc::Int = 2

    """
    Charged species list. Default: [1...nc]
    """
    cspecies::Vector{Int}=collect(1:nc)
    
    "Number of surface species"
    na::Int = 0

    """
    Surface species list. Default: [(nc + 1)...(nc + na)]
    """
    sspecies::Vector{Int}=collect((nc + 1):(nc + na))

    "Index of electrostatic potential ``ϕ`` in species list."
    iϕ::Int = maximum(cspecies) + na + 1

    "Index of pressure `p` in species list"
    ip::Int = maximum(cspecies) + na + 2

    "Mobility coefficients ``D_i\\; (i=1…N)``"
    D::Vector{Float64} = fill(2.0e-9 * ufac"m^2/s", maximum(cspecies))

    "Charge numbers of ions ``z_i\\; (i=1…N)``"
    z::Vector{Int} = [(-1)^(i - 1) for i in 1:maximum(cspecies)]

    "Molar weight of solvent ``M_0``"
    M0::Float64 = 18.0153 * ufac"g/mol"

    "Molar weights of ions ``M_i\\; (i=1…N)``"
    M::Vector{Float64} = fill(M0, maximum(cspecies))

    "Molar volume of solvent ``v_0``"
    v0::Float64 = 1 / (55.4 * ufac"M")

    "Molar volumes of ions ``v_i\\; (i=1…N)``"
    v::Vector{Float64} = fill(v0, maximum(cspecies))

    "Solvation numbers of ions ``κ_i\\; (i=1…N)``"
    κ::Vector{Float64} = fill(10.0, maximum(cspecies))

    """
        actcoeff!(γ, c, p, ::ElectrolyteData)

    Activity coefficient function. Write activity coefficients  ``γ_i\\; (i=1…N)`` into vector `γ`.
    Default: [`DGML_gamma!`](@ref). 
    """
    actcoeff!::Tγ = DGML_gamma!

    "Bulk ion concentrations ``c_i^b\\; (i=1…N)`` "
    c_bulk::Vector{Float64} = fill(0.1 * ufac"M", maximum(cspecies))

    "Working electrode boundary number"
    Γ_we::Int = 1

    "Bulk boundary number"
    Γ_bulk::Int = 2

    "Bulk voltage ``ϕ^b``"
    ϕ_bulk::Float64 = 0.0 * ufac"V"

    "Bulk pressure ``p^b``"
    p_bulk::Float64 = 0.0 * ufac"Pa"

    "Temperature ``T``"
    T::Float64 = (273.15 + 25) * ufac"K"

    "Dielectric permittivity of solvent ``ε``"
    ε::Float64 = 78.49

    "Regularized exponential, default: `exp` (unregularized)"
    rexp::Texp = exp

    "Regularized logarithm, default: `log` (unregularized)"
    rlog::Tlog = log

    "Pressure scaling factor. Default: 1.0e9"
    pscale::Float64 = 1.0e9

    "Local electroneutrality switch. Default: false"
    eneutral::Bool = false

    """
    Upwind flux caculation method for ionic species.
    This allows to choose between
    -  [`μex_flux!`](@ref) (default, strongly preferrable): excess chemical potential (SEDAN) scheme
    -  [`act_flux!`](@ref): scheme based on reciprocal activity coefficients
    -  [`cent_flux!`](@ref): central scheme
    """
    upwindflux!::Tflux = μex_flux!

    """
    Species weights for norms in solver control.
    """
    weights::Vector{Float64} = [v..., zeros(na)..., 1.0, 0.0]

    """
    Solve for pressure. 

    This is `true` by default. Setting this to `false` can serve two purposes:
    - Use the pressure from the solution of a flow equation
    - Ignore the pressure contribution to the excess chemical potential
    """
    solvepressure::Bool = true

    "Faraday constant ``F`` (fixed)"
    F::Float64 = ph"N_A" * ph"e"

    "Dielectric permittivity of vacuum ``ε_0`` (fixed)"
    ε_0::Float64 = ph"ε_0"

    "Molar gas constant scaled with temperature ``RT`` (derived)"
    RT::Float64 = ph"R" * T

    "Solvated molar mass ratio ``m_i = \\frac{M_i}{M_0} + κ_i \\; (i=1… N)`` (derived)"
    Mrel::Vector{Float64} = M / M0 + κ

    "Solvated molar volume ratio ``\\hat v_i =  \\frac{v_i}{v_0} + κ_i \\; (i=1… N)`` (derived)"
    vrel::Vector{Float64} = v / v0 + κ

    "Solvated molar volume ``\\bar v_i = v_i + κ_i v_0 \\; (i=1… N)`` (derived)"
    barv::Vector{Float64} = v + κ * v0

    "Pressure relevant volume ``\\tilde v_i = \\bar v_i + m_i v_0 \\; (i=1… N)``  (derived)"
    tildev::Vector{Float64} = barv - Mrel * v0

    "Activity coefficients at bulk interface ``γ_i^b= γ_i(c_1^b … c_N^b, p^b) \\; (i=1… N)`` (derived)"
    γ_bulk::Vector{Float64} = ones(maximum(cspecies))

    "Cache for activity coefficient calculation (reserved)"
    γk_cache::Tcache = DiffCache(zeros(maximum(cspecies)), 10 * maximum(cspecies))

    "Cache for activity coefficient calculation (reserved)"
    γl_cache::Tcache = DiffCache(zeros(maximum(cspecies)), 10 * maximum(cspecies))

    """
    Working electrode voltage ``ϕ_{we}`` (reserved)
    Used by sweep algorithms to pass boundary data.
    """
    ϕ_we::Float64 = 0.0 * ufac"V"

    """
    Edge velocity projection (reserved).
    """
    edgevelocity::Union{Float64, Vector{Float64}} = 0.0

    """
    Scheme parameter. Deprecated and disabled. 
    Use `upwindflux!` instead to change the discretization scheme.
    """
    scheme::Symbol = :deprecated
end

function Base.setproperty!(this::ElectrolyteData, key::Symbol, value)
    if key == :scheme
        @warn """ Setting ElectrolyteData.scheme  is deprecated and has been disabled. 
            Use `upwindflux!` instead to change the discretization scheme.
        """
        return nothing
    elseif key == :edgevelocity
        return Base.setfield!(this, key, value)
    else
        return Base.setfield!(this, key, convert(typeof(getfield(this, key)), value))
    end
end

function Base.show(io::IOContext{Base.TTY}, this::ElectrolyteData)
    return showstruct(io, this)
end

function Base.show(io::IOContext{Base.IOBuffer}, this::ElectrolyteData)
    return showstruct(io, this)
end

"""
    update_derived!(electrolyte::ElectrolyteData)

Update derived electrolyte data. This needs to be called in order to update fields
`Mrel`, `vrel`, `barv`, `tildev`, `RT`, `γ_bulk` of `ElectrolyteData` after changing some
of the other parameters.

Called on the passed electrolyte data by [`PNPSystem`](@ref),  [`PBSystem`](@ref),  [`dlcapsweep`](@ref),  [`ivsweep`](@ref),  [`cvsweep`](@ref).
"""
function update_derived!(electrolyte::ElectrolyteData)
    (; M, M0, κ, T, v, v0, T, c_bulk, p_bulk, actcoeff!) = electrolyte
    electrolyte.Mrel .= M / M0 + κ
    electrolyte.vrel .= v / v0 + κ
    electrolyte.barv .= v + κ * v0
    electrolyte.tildev .= electrolyte.barv - electrolyte.Mrel * v0
    electrolyte.RT = ph"R" * T
    actcoeff!(electrolyte.γ_bulk, c_bulk, p_bulk, electrolyte)
    return electrolyte
end

function update_derived!(electrolytes::Vector{Td}) where {Td <: AbstractElectrolyteData}
    map(update_derived!, electrolytes)
end

"""
    solvepressure(electrolyte)

Check if pressure is to be solved for using the momentum balance equation, or if the pressure is 
obtained from the (Navier)-Stokes solver.
"""
solvepressure(electrolyte) = electrolyte.solvepressure

"""
    value_or_getindex(v,i)

Return number or i-th component of `Union{Number, Vector}`.
"""
value_or_getindex(v::Number, i) = v
value_or_getindex(v::Vector, i) = v[i]

"""
    edgevelocity(electrolyte, iedge)

Obtain the velocity projection onto a simplex edge of index i (normals to Voronoi cell boundaries).
"""
function edgevelocity(electrolyte, i)
    return value_or_getindex(electrolyte.edgevelocity, i)
end


"""
    dlcap0(electrolyte)

Return double layer capacitance at zero voltage for symmetric binary electrolyte.
The molarity is defined by the bulk concentration value ``c_1^b``.
```math
    c_{dl}=2 \\frac{εε_0F^2c_1^b}{RT}
```

### Example
```jldoctest
using LessUnitful
ely = ElectrolyteData(c_bulk=fill(0.01ufac"mol/dm^3",2))
round(dlcap0(ely),sigdigits=5) |> u"μF/cm^2"
# output

22.847 μF cm^-2
```
"""
function dlcap0(data::AbstractElectrolyteData)
    return sqrt(2 * data.ε * data.ε_0 * data.F^2 * data.c_bulk[1] / (data.RT))
end

"""
    debyelength(electrolyte)

Return debye length for symmetric binary electrolyte. The molarity is defined
by the bulk concentration value ``c_1^b``.
```math
    l_{debye}=\\sqrt{RT\\frac{εε_0}{F^2c_1^b}}
```

```jldoctest
using LessUnitful
ely = ElectrolyteData(c_bulk=fill(0.01ufac"mol/dm^3",2))
round(debyelength(ely),sigdigits=5) |> u"nm"
# output

4.3018 nm
```
"""
debyelength(data) = sqrt(data.ε * data.ε_0 * data.RT / (data.F^2 * data.c_bulk[1]))

"""
    chargedensity(c,electrolyte)

Calculate charge density from vector of concentrations  (in one grid point).
"""
function chargedensity(u::AbstractVector, electrolyte::AbstractElectrolyteData)
    q = zero(eltype(u))
    for ic in electrolyte.cspecies
        q += u[ic] * electrolyte.z[ic]
    end
    return q * electrolyte.F
end

"""
    chargedensity!(q, sol,electrolyte)

Calculate charge density from solution (on the whole grid), writing the result into `q` and returning this vector.
"""
function chargedensity!(q::AbstractVector, u::AbstractMatrix, electrolyte::AbstractElectrolyteData)
    nnodes = size(u, 2)
    for i in 1:nnodes
        @views q[i] = chargedensity(u[:, i], electrolyte)
    end
    return q
end

"""
    chargedensity(sol,electrolyte)

Calculate charge density vector from solution (on whole grid)
"""
function chargedensity(u::AbstractMatrix, electrolyte::AbstractElectrolyteData)
    return chargedensity!(zeros(size(u, 2)), u, electrolyte)
end


"""
    chargedensity(tsol,electrolyte)

Calculate time dependent charge densities from  from time/voltage dependent solution solution (on whole grid).
Returns a [`VoronoiFVM.TransientSolution`](https://wias-pdelib.github.io/VoronoiFVM.jl/stable/solutions/#VoronoiFVM.TransientSolution).
"""
function chargedensity(tsol::TransientSolution, electrolyte)
    nv = length(tsol.t)
    nx = size(tsol.u[1], 2)
    charges = [ reshape(chargedensity(tsol.u[i], electrolyte), (1, nx)) for i in 1:nv]
    return TransientSolution(charges, tsol.t)
end


"""
	c0_barc(u,electrolyte)

Calculate solvent concentration ``c_0`` and summary concentration ``\\bar c`` from vector of concentrations `c`
using the incompressibility constraint (assuming ``κ_0=0``):
```math
 \\sum_{i=0}^N c_i (v_i + κ_iv_0) =1
```

This gives

```math
 c_0v_0=1-\\sum_{i=1}^N c_i (v_i+ κ_iv_0)
```

```math
c_0= 1/v_0 - \\sum_{i=1}^N c_i(\\frac{v_i}{v_0}+κ_i)
```

Then we can calculate 
```math
 \\bar c= \\sum_{i=0}^N c_i
```
"""
function c0_barc(c, electrolyte)
    (; v0, vrel, cspecies) = electrolyte
    T = eltype(c)
    c0 = one(T) / v0
    barc = zero(T)
    for ic in cspecies
        barc += c[ic]
        c0 -= c[ic] * vrel[ic]
    end
    barc += c0
    return c0, barc
end


"""
       solventconcentration(U::Array, electrolyte)

Calculate vector of solvent concentrations from solution array.
"""
function solventconcentration(U::Array, electrolyte)
    (; v0, vrel, nc) = electrolyte
    @views c0 = similar(U[1, :])
    c0 .= 1.0 / v0
    for ic in 1:(nc)
        @views  c0 .-= U[ic, :] .* vrel[ic]
    end
    return c0
end

"""
     chemical_potential(c, barc, p, barv, electrolyte)

Calculate chemical potential of species with concentration c

```math
        μ = \\bar v(p-p_{ref}) + RT\\log \\frac{c}{\\bar c}
```
"""
function chemical_potential(c, barc, p, barv, electrolyte)
    (; rlog, RT, pscale, p_bulk) = electrolyte
    return rlog(c / barc) * RT + barv * (pscale * p - p_bulk)
end

"""
    chemical_potentials!(μ,u,electrolyte)

Calculate chemical potentials from concentrations.

Input:
  -  `μ`: memory for result (will be filled)
  -  `u`: local solution vector (concentrations, potential, pressure)
Returns `μ0, μ`: chemical potential of solvent and chemical potentials of ions.


```jldoctest
using LessUnitful
ely = ElectrolyteData(c_bulk = fill(0.01ufac"mol/dm^3", 2))
μ0, μ = chemical_potentials!([0.0, 0.0], vcat(ely.c_bulk, [0, 0]), ely)
round(μ0, sigdigits = 5), round.(μ, sigdigits = 5)
# output

(-0.89834, [-21359.0, -21359.0])
```

!!! note "Breaking change between v0.2 and 0.3"
    This function has been corrected in version 0.3.  Before it used the molar
    volume ``v`` instead of the effective molar volume ``\\bar v = v + κv_0``
    to calculate the ion chemical potentials. Examples with κ=0 are not affected.
"""
function chemical_potentials!(μ, u::AbstractVector, electrolyte::AbstractElectrolyteData)
    (; ip, v0, v, cspecies, v, v0, barv) = electrolyte
    c0, barc = c0_barc(u, electrolyte)
    μ0 = chemical_potential(c0, barc, u[ip], v0, electrolyte)
    for i in cspecies
        μ[i] = chemical_potential(u[i], barc, u[ip], barv[i], electrolyte)
    end
    return μ0, μ
end

"""
    chemical_potentials(solution, electrolyte)

Calculate chemical potentials from solution.
Returns `μ0, μ`,  where `μ0` is the vector of solvent chemical potentials,
and `μ` is the `nc×nnodes` matrix of solute chemical potentials.
"""
function chemical_potentials(u::AbstractMatrix, data::AbstractElectrolyteData)
    nn = size(u, 2)
    μ = zeros(maximum(data.cspecies), nn)
    μ0 = zeros(nn)
    for i in 1:nn
        μ0[i], _ = chemical_potentials!(view(μ, :, i), u[:, i], data)
    end
    return μ0, μ
end

"""
    electrochemical_potentials(solution, electrolyte)

Calculate electrochemical potentials from solution.
and return `nc×nnodes` matrix of solute electrochemical potentials measured in Volts.
"""
function electrochemical_potentials(u::AbstractMatrix, data::AbstractElectrolyteData)
    μ0, μ = chemical_potentials(u, data)
    (; iϕ, cspecies, z, F, RT, M, M0) = data
    nn = size(u, 2)
    for i in 1:nn
        for ic in cspecies
            μ[ic, i] /= F
            μ[ic, i] += z[ic] * u[iϕ, i] - μ0[i] * M[ic] / (M0 * F)
        end
    end
    return μ
end


"""
    rrate(R0,β,A)
    rrate(R0,β,A, electrolyte)

Reaction rate expression

    rrate(R0,β,A)=R0*(exp(-β*A) - exp((1-β)*A))
"""
function rrate(R0, β, A)
    return R0 * (exp(-β * A) - exp((1 - β) * A))
end

function rrate(R0, β, A, electrolyte)
    (; rexp) = electrolyte
    return R0 * (rexp(-β * A) - rexp((1 - β) * A))
end


"""
    wnorm(u,w,p)

Weighted norm with respect to columns
"""
function wnorm(u, w, p)
    @views norms = [w[i] * LinearAlgebra.norm(u[i, :], p) for i in 1:size(u, 1)]
    return LinearAlgebra.norm(norms, p)
end


"""
    isincompressible(cx::Vector,electrolyte)

Check if the concentration vector fulfills incompressibility constraint, including
the fact that the solvent concentration is nonnegative
"""
function isincompressible(cx::Vector, electrolyte)
    (; v0, v, κ) = electrolyte
    c0, barc = c0_barc(cx, electrolyte)
    cv = c0 * v0
    for i in electrolyte.cspecies
        cv += cx[i] * (v[i] + κ[i] * v0)
    end
    return c0 ≥ 0 && cv ≈ 1.0
end

"""
    isincompressible(tsol::TransientSolution,electrolyte)

Check for incompressibility of transient solution
"""
function isincompressible(tsol::TransientSolution, electrolyte)
    return all(cx -> isincompressible(cx, electrolyte), [u[:, i] for u in tsol.u, i in size(tsol, 2)])
end

"""
    iselectroneutral(cx::Vector,electrolyte)

Check for electroneutrality of concentration vector
"""
function iselectroneutral(cx, electrolyte)
    return isapprox(cx[electrolyte.cspecies]' * electrolyte.z[electrolyte.cspecies], 0; atol = 1.0e-12)
end

"""
    iselectroneutral(tsol::TransientSolution,electrolyte)

Check for electroneutrality of transient solution
"""
function iselectroneutral(tsol::TransientSolution, electrolyte)
    return all(cx -> iselectroneutral(cx, electrolyte), [u[:, i] for u in tsol.u, i in size(tsol, 2)])
end


"""
    bulkbcondition(f,u,bnode,electrolyte; region = data.Γ_bulk)

Bulk boundary condition for electrolyte: set potential, pressure and concentrations to bulk values.
"""
function bulkbcondition(f, u, bnode, electrolyte; region = electrolyte.Γ_bulk)
    (; iϕ, ip, cspecies, ϕ_bulk, p_bulk, c_bulk) = electrolyte
    if bnode.region == region
        boundary_dirichlet!(f, u, bnode; species = iϕ, region, value = ϕ_bulk)
        boundary_dirichlet!(f, u, bnode; species = ip, region, value = p_bulk)
        for ic in cspecies
            boundary_dirichlet!(f, u, bnode; species = ic, region, value = c_bulk[ic])
        end
    end
    return nothing
end
