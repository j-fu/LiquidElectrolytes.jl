"""
    pbreaction(f, u, node, electrolyte)

Reaction expression for Poisson-Boltzmann
"""
function pbreaction(f, u, node, electrolyte)
    (; ip, iϕ, v0, v, M0, z, M, κ, F, RT, nc, pscale, p_bulk, c_bulk,
     γ!,γk_cache, γl_cache
     ) = electrolyte
    p = u[ip] * pscale - p_bulk
    ϕ = u[iϕ]
    γ=get_tmp(γk_cache, u)
    γ_bulk=get_tmp(γl_cache, u)
    γ!(γ_bulk, c_bulk, p_bulk, electrolyte)
    γ!(γ, u, p, electrolyte)
    q = zero(eltype(u))
    for ic in 1:nc
        f[ic] = u[ic] - c_bulk[ic] * rexp(-z[ic] * F * ϕ / RT)*γ_bulk[ic]/γ[ic]
        q += z[ic] * u[ic]
    end
    f[iϕ] = -F * q
    f[ip] = 0
    return
end


"""
    pbflux(f, u, edge, electrolyte)

Flux expression for Poisson-Boltzmann
"""
function pbflux(f, u, edge, electrolyte)
    (; ε_0, ε, iϕ, ip, nc, z, F, pscale) = electrolyte
    dϕ= u[iϕ, 1] - u[iϕ, 2]

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
        celldata = ElectrolyteData(),
        bcondition = (f, u, n, e) -> nothing,
        kwargs...
    )

    sys= VoronoiFVM.System(
        grid;
        data = celldata,
        flux = pbflux,
        reaction = pbreaction,
        bcondition,
        species = [1:(celldata.nc)..., celldata.iϕ, celldata.ip],
        kwargs...
    )

    return PBSystem(sys)
end
