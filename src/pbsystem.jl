"""
    pbreaction(f, u, node, electrolyte)

Reaction expression for Poisson-Boltzmann
"""
function pbreaction(f, u, node, electrolyte)
    (; ip, iϕ, v0, v, M0, z, M, κ, F, RT, nc, pscale, p_bulk, c_bulk) = electrolyte
    c0, bar_c = c0_barc(u, electrolyte)
    c0_bulk, barc_bulk = c0_barc(c_bulk, electrolyte)
    p = u[ip] * pscale - p_bulk
    ϕ = u[iϕ]
    q = zero(eltype(u))
    for ic in 1:nc
        if v[ic]==0
            f[ic] = u[ic] - c_bulk[ic] * rexp(-z[ic] * F * ϕ / RT)
        else
            Mrel = M[ic] / M0 + κ[ic]
            barv = v[ic] + κ[ic] * v0
            tildev = barv - Mrel * v0
            γ_bulk = barc_bulk * (c0_bulk / barc_bulk)^Mrel
            γ = bar_c * rexp( -tildev * p / RT) * (c0 / bar_c)^Mrel
            f[ic] = u[ic] - c_bulk[ic] * rexp(-z[ic] * F * ϕ / RT) * γ / γ_bulk
        end
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
    f[iϕ] = ε * ε_0 * (u[iϕ, 1] - u[iϕ, 2])
    qk = zero(eltype(u))
    ql = zero(eltype(u))
    for ic in 1:nc
        qk += z[ic] * u[ic, 1]
        ql += z[ic] * u[ic, 2]
    end
    f[ip] = (u[ip, 1] - u[ip, 2]) + F * (u[iϕ, 1] - u[iϕ, 2]) * (qk + ql) / 2 * pscale
    return
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
        grid;
        celldata = ElectrolyteData(),
        bcondition = (f, u, n, e) -> nothing,
        kwargs...
    )

    return VoronoiFVM.System(
        grid;
        data = celldata,
        flux = pbflux,
        reaction = pbreaction,
        bcondition,
        species = [1:(celldata.nc)..., celldata.iϕ, celldata.ip],
        kwargs...
    )
end
