"""
    pbspacecharge(φ, p, electorlyte)

Space charge expression for Poisson-Boltzmann
"""
function pbspacecharge(φ, p, data)
    c0_bulk, barc_bulk = c0_barc(data.c_bulk, data)
    pscaled = (p * data.pscale - data.p_bulk)
    c0 = c0_bulk * exp(-data.v0 * pscaled / (data.RT))
    sumyz = zero(p)
    sumyv = data.v0 * c0
    for α in 1:data.nc
        barv = data.v[α] + data.κ[α] * data.v0
        η_p = barv * pscaled
        η_φ = data.z[α] * data.F * (φ - data.ϕ_bulk)
        y = data.c_bulk[α] * exp(-(η_φ + η_p) / (data.RT))
        sumyz += data.z[α] * y
        sumyv += barv * y
    end
    return data.F * sumyz / sumyv
end

"""
    pbreaction(f, u, node, electrolyte)

Reaction expression for Poisson-Boltzmann
"""
function pbreaction(f, u, node, electrolyte)
    iϕ, ip = 1, 2
    ## Charge density
    f[iϕ] = -pbspacecharge(u[iϕ], u[ip], electrolyte)
    f[ip] = 0
    return 
end

"""
    pbflux(f, u, edge, electrolyte)

Flux expression for Poisson-Boltzmann
"""
function pbflux(f, u, edge, data)
    iϕ, ip = 1, 2
    f[iφ] = data.ε * data.ε_0 * (u[iφ, 1] - u[iφ, 2])
    if iszero(data.v)
        qavg = 0
    else
        q1 = pbspacecharge(u[iφ, 1], u[ip, 1], data)
        q2 = pbspacecharge(u[iφ, 2], u[ip, 2], data)
        qavg = (q1 + q2) / 2
    end
    f[ip] = (u[ip, 1] - u[ip, 2]) + (u[iφ, 1] - u[iφ, 2]) * qavg / data.pscale 
    return 
end


"""
    PBSystem(grid;
             celldata=ElectrolyteData(),
             bcondition=default_bcondition,
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
        bcondition = default_bcondition,
        kwargs...
    )

    return sys = VoronoiFVM.System(
        grid;
        data = celldata,
        flux = pbflux,
        reaction = pbreaction,
        bcondition,
        species = [1, 2],
        kwargs...
    )
end
