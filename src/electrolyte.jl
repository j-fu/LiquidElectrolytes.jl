"""
$(TYPEDEF)

Abstract super type for electrolytes
"""
abstract type AbstractElectrolyteData{Q} end

"""
$(TYPEDEF)

Data for electrolyte. It is defined using `Base.@kwdef`
has keyword constructors like
```julia
    ElectrolyteData(nc=3,z=[-1,2,1])
```

$(TYPEDFIELDS)
"""
@kwdef mutable struct ElectrolyteData{Q} <: AbstractElectrolyteData{Q}
    "Number of ionic species."
    nc::Int = 2

    "Number of surface species"
    na::Int = 0

    "Potential index in species list."
    iϕ::Int = nc + na + 1

    "Pressure index in species list"
    ip::Int = nc + na + 2

    "Mobility coefficient"
    D::Vector{Q} = fill(2.0e-9 * u"m^2/s", nc) |> Vector

    "Charge numbers of ions"
    z::Vector{Int} = [(-1)^(i - 1) for i = 1:nc]

    "Molar weight of solvent"
    M0::Q = 18.0153 * u"g/mol"

    "Molar weight of ions"
    M::Vector{Q} = fill(M0, nc)  |> Vector

    "Molar volume of solvent"
    v0::Q = 1 / (55.4 * u"mol/dm^3")

    "Molar volumes of ions"
    v::Vector{Q} = fill(v0, nc)  |> Vector

    "Solvation numbers"
    κ::Vector{Float64} = fill(10.0, nc)  |> Vector

    "Bulk ion concentrations"
    c_bulk::Vector{Q} = fill(0.1 * u"mol/L", nc)  |> Vector

    "Bulk voltage"
    ϕ_bulk::Q = 0.0 * u"V"

    "Bulk pressure"
    p_bulk::Q = 0.0 * u"Pa"

    "Bulk boundary number"
    Γ_bulk::Int = 2

    "Working electrode voltage"
    ϕ_we::Q = 0.0 * u"V"

    "Working electrode  boundary number"
    Γ_we::Int = 1

    "Temperature"
    T::Q = (273.15 + 25) * u"K"

    "Molar gas constant scaled with temperature"
    RT::Q = Constants.R * T

    "Faraday constant"
    F::Q = Constants.N_A*Constants.e

    "Dielectric permittivity of solvent"
    ε::Float64 = 78.49

    "Dielectric permittivity of vacuum"
    ε_0::Q = Constants.eps_0

    "Pressure scaling factor"
    pscale::Q = 1.0e9*u"Pa"

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
    Regularization parameter used in [`rlog`](@ref)
    """
    epsreg::Float64 = 1.0e-20
    """
    Species weights for norms in solver control.
    """
    weights::Vector{Float64} = [ones(nc)..., zeros(na)..., 1.0, 0.0]
end
ustrip(q::Number)=q
ustrip(s::Symbol)=s
ustrip(vq::AbstractVector)=[ustrip(q) for q in vq]
ustrip(this::AbstractElectrolyteData{Q}) where {Q} = ElectrolyteData{Float64}([ustrip(getfield(this,name)) for name in fieldnames(ElectrolyteData)]...)


function Base.show(io::IO, this::ElectrolyteData)
    showstruct(io, this)
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
    sqrt(2 * data.ε * data.ε_0 * data.F^2 * data.c_bulk[1] / (data.RT))
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
    charge(c,electrolyte)

Calculate charge from vector of concentrations
"""
function charge(u, electrolyte::AbstractElectrolyteData)
    q = zero(eltype(u))
    for ic = 1:(electrolyte.nc)
        q += u[ic] * electrolyte.z[ic]
    end
    q * electrolyte.F
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
    for ic = 1:(electrolyte.nc)
        barc += c[ic]
        c0 -= c[ic] * vrel(ic, electrolyte)
    end
    barc += c0
    c0, barc
end

"""
    rlog(u, electrolyte)

Calls rlog(u;eps=electrolyte.epsreg)
"""
rlog(x, electrolyte::AbstractElectrolyteData) = rlog(x; eps = electrolyte.epsreg)

"""
    rlog(u; eps=1.0e-20)

Regularized logarithm. Smooth linear continuation for `x<eps`.
This means we can calculate a "logarithm"  of a small negative number.
"""
function rlog(x; eps = 1.0e-20)
    if x < eps
        return log(eps) + (x - eps) / eps
    else
        return log(x)
    end
end

"""
       solventconcentration(U::Array, electrolyte)

Calculate vector of solvent concentrations from solution array.
"""
function solventconcentration(U::Array, electrolyte)
    c0 = similar(U[1, :])
    c0 .= 1.0 / electrolyte.v0
    for ic = 1:(electrolyte.nc)
        c0 -= U[ic, :] .* vrel(ic, electrolyte)
    end
    c0
end

@doc raw"""
        chemical_potential(c, barc, p, v, electrolye)

Calculate chemical potential of species with concentration c

```math
        μ = v(p-p_{ref}) + RT\log \frac{c}{\bar c}
```
"""
chemical_potential(c, barc, p, v, data) = rlog(c / barc, data) * data.RT + v * data.pscale * (p - data.p_bulk)

"""
    chemical_potentials!(μ,u,electrolyte)

Calculate chemical potentials from concentrations.

Input:
  -  `μ`: memory for result (will be filled)
  -  `u`: local solution vector (concentrations, potential, pressure)
Returns `μ0, μ`: chemical potential of solvent and chemical potentials of ions.


```jldoctest
using LessUnitful
ely = ElectrolyteData(c_bulk=fill(0.01ufac"mol/dm^3",2))
μ0,μ=chemical_potentials!([0.0,0.0],vcat(ely.c_bulk,[0,0]),ely)
round(μ0,sigdigits=5),round.(μ,sigdigits=5)
# output

(-0.89834, [-21359.0, -21359.0])
```
"""
function chemical_potentials!(μ, u, data::AbstractElectrolyteData)
    (; ip, pscale, RT, v0, v, nc) = data
    c0, barc = c0_barc(u, data)
    μ0 = chemical_potential(c0, barc, u[ip], data.v0, data)
    for i = 1:nc
        μ[i] = chemical_potential(u[i], barc, u[ip], data.v[i], data)
    end
    μ0, μ
end

"""
    rexp(x;trunc=500.0)

Regularized exponential. Linear continuation for `x>trunc`,  
returns 1/rexp(-x) for `x<-trunc`.
"""
function rexp(x; trunc = 20.0)
    if x < -trunc
        1.0 / rexp(-x; trunc)
    elseif x <= trunc
        exp(x)
    else
        exp(trunc) * (x - trunc + 1)
    end
end

"""
    rrate(R0,β,A)

Reaction rate expression

    rrate(R0,β,A)=R0*(exp(-β*A) - exp((1-β)*A))
"""
rrate(R0, β, A) = R0 * (rexp(-β * A) - rexp((1 - β) * A))

"""
    wnorm(u,w,p)

Weighted norm with respect to columns
"""
function wnorm(u, w, p)
    @views norms = [w[i] * LinearAlgebra.norm(u[i, :], p) for i = 1:size(u, 1)]
    LinearAlgebra.norm(norms, p)
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
    for i = 1:(celldata.nc)
        cv += cx[i] * (v[i] + κx[i] * v0)
    end
    cv ≈ 1.0
end

"""
    isincompressible(tsol::TransientSolution,celldata)

Check for incompressibility of transient solution
"""
function isincompressible(tsol::TransientSolution, celldata)
    all(cx -> isincompressible(cx, celldata), [u[:, i] for u in tsol.u, i in size(tsol, 2)])
end

"""
    iselectroneutral(cx::Vector,celldata)

Check for electroneutrality of concentration vector
"""
function iselectroneutral(cx, celldata)
    isapprox(cx[1:(celldata.nc)]' * celldata.z, 0; atol = 1.0e-12)
end

"""
    iselectroneutral(tsol::TransientSolution,celldata)

Check for electorneutrality of transient solution
"""
function iselectroneutral(tsol::TransientSolution, celldata)
    all(cx -> iselectroneutral(cx, celldata), [u[:, i] for u in tsol.u, i in size(tsol, 2)])
end
