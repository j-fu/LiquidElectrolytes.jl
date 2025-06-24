"""l
    pbreaction!(f, u, node, electrolyte)

Reaction expression for Poisson-Boltzmann
"""
function pbreaction!(f, u, node, electrolyte)
    (;
        ip, iϕ, cspecies,
        z, F, RT, pscale,
        c_bulk, ϕ_bulk, p_bulk,
        actcoeff!, γk_cache, rexp, γ_bulk,
    ) = electrolyte
    p = u[ip] * pscale - p_bulk
    ϕ = u[iϕ]
    γ = get_tmp(γk_cache, u)
    actcoeff!(γ, u, p, electrolyte)
    q = zero(eltype(u))
    for ic in cspecies
        f[ic] = u[ic] * γ[ic] - c_bulk[ic] * rexp(z[ic] * F * (ϕ_bulk - ϕ) / RT) * γ_bulk[ic]
        q += z[ic] * u[ic]
    end
    f[iϕ] = -F * q
    return
end


"""
    pbflux!(f, u, edge, electrolyte)

Flux expression for Poisson-Boltzmann
"""
function pbflux!(f, u, edge, electrolyte)
    (; ε_0, ε, iϕ, ip, pscale) = electrolyte
    dϕ = u[iϕ, 1] - u[iϕ, 2]
    f[iϕ] = ε * ε_0 * dϕ
    @views qk, ql = chargedensity(u[:, 1], electrolyte), chargedensity(u[:, 2], electrolyte)
    f[ip] = u[ip, 1] - u[ip, 2] + (qk + ql) * dϕ / (2 * pscale)
    return
end

struct PBSystem <: AbstractElectrochemicalSystem
    vfvmsys::VoronoiFVM.System
end

"""
    PBSystem(grid;
                 celldata=ElectrolyteData(),
             bcondition=(f, u, n, e)-> nothing,
             kwargs...)

Create VoronoiFVM system generalized Poisson-Boltzmann. Input:
- `grid`: discretization grid
- `celldata`: composite struct containing electrolyte data
- `bcondition`: boundary condition
- `kwargs`: Keyword arguments of VoronoiFVM.System
"""
function PBSystem(
        grid::ExtendableGrid;
        celldata::ElectrolyteData = ElectrolyteData(),
        bcondition = (f, u, n, e) -> nothing,
        kwargs...
    )
    update_derived!(celldata)
    sys = VoronoiFVM.System(
        grid;
        data = celldata,
        flux = pbflux!,
        reaction = pbreaction!,
        bcondition,
        species = union(celldata.cspecies, [celldata.ip, celldata.iϕ]),
        kwargs...
    )

    return PBSystem(sys)
end
