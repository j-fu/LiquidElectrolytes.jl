"""
   Structures & methods from v1.0 still used in IonChannelProject
"""
abstract type AbstractLegacyElectrolyteData <: AbstractElectrolyteData end


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
@kwdef mutable struct LegacyElectrolyteData <: AbstractLegacyElectrolyteData
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

end

"""
    sflux(ic,dϕ,ck,cl,γk,γl,bar_ck,bar_cl,electrolyte; evelo=0.0)

 Sedan flux,  see Gaudeul/Fuhrmann 2022

Appearantly first described by Yu, Zhiping  and Dutton, Robert, SEDAN III, www-tcad.stanford.edu/tcad/programs/sedan3.html

 see also the 198? Fortran code available via
 https://web.archive.org/web/20210518233152/http://www-tcad.stanford.edu/tcad/programs/oldftpable.html

Verification calculation is in the paper.
"""
function sflux(ic, dϕ, ck, cl, γk, γl, bar_ck, bar_cl, electrolyte; evelo = 0.0)
    (; D, z, F, RT) = electrolyte
    bp, bm = fbernoulli_pm(z[ic] * dϕ * F / RT + dμex(γk, γl, electrolyte) / RT - evelo / D[ic])
    return D[ic] * (bm * ck - bp * cl)
end

#=

B(b-a)/B(a-b)= (b-a)*(exp(a-b)-1)/ (a-b)*(exp(b-a)-1)
             = (1-exp(a-b))/(exp(b-a)-1)
             = exp(-b)*(exp(b)-exp(a)) / exp(-a)*(exp(b)-exp(a))
             = exp(a)/exp(b)

c/barc= C
ck/cl = bp/bm = exp(z ϕk*F/RT + μex_k/RT)/exp(z ϕl*F/RT + μex_l/RT)
=#

"""
    aflux(ic,dϕ,ck,cl,γk,γl,bar_ck,bar_cl,electrolyte; evelo=0)

Flux expression based on  activities, see Fuhrmann, CPC 2015
??? Do we need to divide the velocity by the inverse activity coefficient ?

"""
function aflux(ic, dϕ, ck, cl, γk, γl, bar_ck, bar_cl, electrolyte; evelo = 0.0)
    (; D, z, F, RT) = electrolyte
    Dx = D[ic] * (1 / γk + 1 / γl) / 2
    bp, bm = fbernoulli_pm(z[ic] * dϕ * F / RT - evelo / D[ic])
    return Dx * (bm * ck * γk - bp * cl * γl)
end

#=
ck/cl= bp/betaK  / bm/betal
     =  exp(z\phik *F/RT + \muexK/RT)/ exp(z\phil *F/RT + \muexL/RT)
=#

"""
    cflux(ic,dϕ,ck,cl,γk,γl,bar_ck,bar_cl,electrolyte; evelo = 0)

Flux expression based on central differences, see Gaudeul/Fuhrmann 2022
"""
function cflux(ic, dϕ, ck, cl, γk, γl, bar_ck, bar_cl, electrolyte; evelo = 0.0)
    (; D, z, F, RT) = electrolyte
    μk = rlog(ck) * RT
    μl = rlog(cl) * RT
    return D[ic] * 0.5 * (ck + cl) * ((μk - μl + dμex(γk, γl, electrolyte) + z[ic] * F * dϕ) / RT - evelo / D[ic])
end

function pnpstorage(f, u, node, electrolyte::AbstractLegacyElectrolyteData)
    f[electrolyte.iϕ] = zero(eltype(u))
    f[electrolyte.ip] = zero(eltype(u))
    for ic in 1:(electrolyte.nc)
        f[ic] = u[ic]
    end
    return
end

@inline function dμex(γk, γl, electrolyte::AbstractLegacyElectrolyteData)
    return (rlog(γk) - rlog(γl)) * (electrolyte.RT)
end

rlog(x::Any) = Base.log(x)

@doc raw"""
	vrel(ic,electrolyte)

Calculate relative (wrt. solvent) molar volume of i-th species ``v_{i,rel}=κ_i+\frac{v_i}{v_0}``.
"""
vrel(ic, electrolyte) = electrolyte.v[ic] / electrolyte.v0 + electrolyte.κ[ic]

"""
	c0_barc(u,electrolyte)

Calculate solvent concentration ``c_0`` and summary concentration ``\\bar c`` from vector of concentrations `c`
using the incompressibility constraint (assuming ``κ_0=0``):
"""
function c0_barc(c, electrolyte::AbstractLegacyElectrolyteData)
    c0 = one(eltype(c)) / electrolyte.v0
    barc = zero(eltype(c))
    for ic in 1:(electrolyte.nc)
        barc += c[ic]
        c0 -= c[ic] * vrel(ic, electrolyte)
    end
    barc += c0
    return c0, barc
end
