"""
$(TYPEDEF)

Abstract super type for electrolytes.
"""
abstract type AbstractElectrolyteData end

"""
$(TYPEDEF)

Data for electrolyte. It is defined using `Base.@kwdef`
Leading  to keyword constructors like
```julia
    ElectrolyteData(nc=3,z=[-1,2,1])
```

Fields (for default values, see below):

$(TYPEDFIELDS)
"""
@kwdef mutable struct ElectrolyteData{Tlog, Texp} <: AbstractElectrolyteData
    "Number of ionic species."
    nc::Int = 2

    "Number of surface species"
    na::Int = 0

    "Potential index in species list."
    iϕ::Int = nc + na + 1

    "Pressure index in species list"
    ip::Int = nc + na + 2

    "Mobility coefficient"
    D::Vector{Float64} = fill(2.0e-9 * ufac"m^2/s", nc)

    "Charge numbers of ions"
    z::Vector{Int} = [(-1)^(i - 1) for i in 1:nc]

    "Molar weight of solvent"
    M0::Float64 = 18.0153 * ufac"g/mol"

    "Molar weights of ions"
    M::Vector{Float64} = fill(M0, nc)

    "Molar volume of solvent"
    v0::Float64 = 1 / (55.4 * ufac"M")

    "Molar volumes of ions"
    v::Vector{Float64} = fill(v0, nc)

    "Solvation numbers of ions"
    κ::Vector{Float64} = fill(10.0, nc)

    "Bulk ion concentrations"
    c_bulk::Vector{Float64} = fill(0.1 * ufac"M", nc)

    "Bulk voltage"
    ϕ_bulk::Float64 = 0.0 * ufac"V"

    "Bulk pressure"
    p_bulk::Float64 = 0.0 * ufac"Pa"

    "Bulk boundary number"
    Γ_bulk::Int = 2

    "Working electrode voltage"
    ϕ_we::Float64 = 0.0 * ufac"V"

    "Working electrode boundary number"
    Γ_we::Int = 1

    "Temperature"
    T::Float64 = (273.15 + 25) * ufac"K"

    "Molar gas constant scaled with temperature"
    RT::Float64 = ph"R" * T

    "Faraday constant"
    F::Float64 = ph"N_A" * ph"e"

    "Dielectric permittivity of solvent"
    ε::Float64 = 78.49

    "Dielectric permittivity of vacuum"
    ε_0::Float64 = ph"ε_0"

    "Pressure scaling factor"
    pscale::Float64 = 1.0e9

    "Local electroneutrality switch"
    eneutral::Bool = false

    """
    [Flux caculation scheme](@id fluxes)
    This allows to choose between
    - `:μex` (default): excess chemical potential (SEDAN) scheme, see [`sflux`](@ref)
    - `:act` : scheme based on reciprocal activity coefficients, see [`aflux`](@ref)
    - `:cent` : central scheme, see [`cflux`](@ref).
    """
    scheme::Symbol = :μex

    """
    Solve for pressure. 

    This is `true` by default. Setting this to `false` can serve two purposes:
    - Take the pressure from the solution of a flow equation
    - Ignore the pressure contribution to the excess chemical potential
    """
    solvepressure::Bool = true

    """
    Species weights for norms in solver control.
    """
    weights::Vector{Float64} = [v..., zeros(na)..., 1.0, 0.0]

    """
    Edge velocity projection.
    """
    edgevelocity::Union{Float64, Vector{Float64}} = 0.0

    """
    Electrolyte model
    """
    model::Symbol = :DGL

    """
    Logarithm function. Default: `Base.log`, but can be replaced by 
    [`RLog`](@ref)
    """
    log::Tlog = Base.log

    """
    Exponential function. Default: `Base.exp`, but can be replaced by 
    [`RExp`](@ref)
    """
    exp::Texp = Base.exp
end

"""
    set_model!(electrolyte, model)

Force the electrolyte data to be consistent to given model. The following
models are supported:
    - `:default`, `:DGL`: Dreyer/Guhlke/Landstorfer model (varying solvation, molar volumes, molar mass, solve for pressure)
    - `:BAO`:  Borukhov/Andelman/Orland (``κ_i=0``, ``M_i=M_0``, `solvepressure=false`)
    - `:DGM`:  Dreyer/Guhlke/Müller (``κ_i=0``, ``v_i=v_0``)
    - `:ρconst`: Constant density (``M_i=M_0v_i/v_0``) (for consistent coupling with stokes)
"""
function set_model!(electrolyte, model)
    if model==:BAO
        electrolyte.κ.=0
        electrolyte.M.=electrolyte.M0
        electrolyte.solvepressure=false
        electrolyte.model=:BAO
    elseif model==:default || model==:DGL
        electrolyte.model=:DGL
    elseif model==:DGM
        electrolyte.κ.=0
        electrolyte.v.=electrolyte.v0
        electrolyte.model=:DGM
    elseif model==:ρconst
        electrolyte.M.=electrolyte.v.*electrolyte.M/electrolyte.v0
        electrolyte.model=:DGM
    else
        error("Unknown electrolyte model: $(model)")
    end
    return nothing
end
    
"""
    solvepressure(electrolyte)

Check if pressure is to be solved for using the pressure Laplace equation derived
from the Navier-Stokes equation in mechanical equilibrium or if the pressure is
obtained from the (Navier)-Stokes solver.
"""
solvepressure(electrolyte) = electrolyte.solvepressure


_evelo(v::Number, i) = v
_evelo(v::Vector, i) = v[i]

"""
    edgevelocity(electrolyte, iedge)

Obtain the velocity projection onto a simplex edge of index i (normals to Voronoi cell boundaries).
"""
function edgevelocity(electrolyte, i)
    return _evelo(electrolyte.edgevelocity, i)
end


function Base.show(io::IO, this::ElectrolyteData)
    return showstruct(io, this)
end

"""
    dlcap0(electrolyte)

Double layer capacitance at zero voltage for symmetric binary electrolyte.

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

Debye length.

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
    for ic in 1:(electrolyte.nc)
        q += u[ic] * electrolyte.z[ic]
    end
    return q * electrolyte.F
end

"""
    chargedensity!(q, sol,electrolyte)

Calculate charge density from solution (on the whole grid), putting the resulti
into `q` and returning this vector.
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

Calculate charge densities from  from time/voltage dependent solution solution (on whole grid)
"""
function chargedensity(tsol::TransientSolution, electrolyte)
    nv = length(tsol.t)
    nx = size(tsol.u[1], 2)
    charges = [ reshape(chargedensity(tsol.u[i], electrolyte), (1, nx)) for i in 1:nv]
    return TransientSolution(charges, tsol.t)
end

@doc raw"""
	vrel(ic,electrolyte)

Calculate relative (wrt. solvent) molar volume of i-th species ``v_{i,rel}=κ_i+\frac{v_i}{v_0}``.
"""
vrel(ic, electrolyte) = electrolyte.v[ic] / electrolyte.v0 + electrolyte.κ[ic]

@doc raw"""
	c0_barc(u,electrolyte)

Calculate solvent concentration ``c_0`` and summary concentration ``\bar c`` from vector of concentrations `c`
using the incompressibility constraint (assuming ``κ_0=0``):
```math
 \sum_{i=0}^N c_i (v_i + κ_iv_0) =1
```

This gives

```math
 c_0v_0=1-\sum_{i=1}^N c_i (v_i+ κ_iv_0)
```

```math
c_0= 1/v_0 - \sum_{i=1}^N c_i(\frac{v_i}{v_0}+κ_i)
```

Then we can calculate 
```math
 \bar c= \sum_{i=0}^N c_i
```
"""
function c0_barc(c, electrolyte)
    c0 = one(eltype(c)) / electrolyte.v0
    barc = zero(eltype(c))
    for ic in 1:(electrolyte.nc)
        barc += c[ic]
        c0 -= c[ic] * vrel(ic, electrolyte)
    end
    barc += c0
    return c0, barc
end



"""
       solventconcentration(U::Array, electrolyte)

Calculate vector of solvent concentrations from solution array.
"""
function solventconcentration(U::Array, electrolyte)
    @views c0 = similar(U[1, :])
    c0 .= 1.0 / electrolyte.v0
    for ic in 1:(electrolyte.nc)
        @views  c0 .-= U[ic, :] .* vrel(ic, electrolyte)
    end
    return c0
end

@doc raw"""
        chemical_potential(c, barc, p, barv, electrolyte)

Calculate chemical potential of species with concentration c

```math
        μ = \bar v(p-p_{ref}) + RT\log \frac{c}{\bar c}
```
"""
chemical_potential(c, barc, p, barv, electrolyte) = electrolyte.log(c / barc) * electrolyte.RT + barv * electrolyte.pscale * (p - electrolyte.p_bulk)

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
function chemical_potentials!(μ, u::AbstractVector, data::AbstractElectrolyteData)
    (; ip, pscale, RT, v0, v, nc) = data
    c0, barc = c0_barc(u, data)
    μ0 = chemical_potential(c0, barc, u[ip], data.v0, data)
    for i in 1:nc
        μ[i] = chemical_potential(u[i], barc, u[ip], data.v[i] + data.κ[i] * data.v0, data)
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
    μ = zeros(data.nc, nn)
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
    (; iϕ, nc, z, F, RT, M, M0) = data
    nn = size(u, 2)
    for i in 1:nn
        for ic in 1:nc
            μ[ic, i] /= F
            μ[ic, i] += z[ic] * u[iϕ, i] - μ0[i] * M[ic] / (M0 * F)
        end
    end
    return μ
end


"""
    rrate(R0,β,A, exp=Base.exp)

Reaction rate expression

    rrate(R0,β,A; exp=Base.exp)=R0*(exp(-β*A) - exp((1-β)*A))
"""
function rrate(R0, β, A, exp::Texp=Base.exp) where  {Texp}
    return R0 * (exp(-β * A) - exp((1 - β) * A))
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
    isincompressible(cx::Vector,celldata)

Check for incompressibility of concentration vector
"""
function isincompressible(cx::Vector, celldata)
    c0, barc = c0_barc(cx, celldata)
    v0 = celldata.v0
    v = celldata.v
    κx = celldata.κ
    cv = c0 * v0
    for i in 1:(celldata.nc)
        cv += cx[i] * (v[i] + κx[i] * v0)
    end
    return cv ≈ 1.0
end

"""
    isincompressible(tsol::TransientSolution,celldata)

Check for incompressibility of transient solution
"""
function isincompressible(tsol::TransientSolution, celldata)
    return all(cx -> isincompressible(cx, celldata), [u[:, i] for u in tsol.u, i in size(tsol, 2)])
end

"""
    iselectroneutral(cx::Vector,celldata)

Check for electroneutrality of concentration vector
"""
function iselectroneutral(cx, celldata)
    return isapprox(cx[1:(celldata.nc)]' * celldata.z, 0; atol = 1.0e-12)
end

"""
    iselectroneutral(tsol::TransientSolution,celldata)

Check for electroneutrality of transient solution
"""
function iselectroneutral(tsol::TransientSolution, celldata)
    return all(cx -> iselectroneutral(cx, celldata), [u[:, i] for u in tsol.u, i in size(tsol, 2)])
end
