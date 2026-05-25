"""
    pnpstorage!(f, u, node, electrolyte)            

Finite volume storage term
"""
function pnpstorage!(f, u, node, electrolyte::AbstractElectrolyteData)
    (; cspecies) = electrolyte
    for ic in cspecies
        f[ic] = u[ic]
    end
    return
end

function pnpstorage!(f, u, node, celldata::AbstractCellData)
    elys = electrolytes(celldata)
    pnpstorage!(f, u, node, elys[node.region])
    return
end

"""
    pnpbstorage!(f, u, node, electrolyte)

Finite volume boundary storage term
"""
function pnpbstorage!(f, u, node, electrolyte)
    (; nc, na, Γ_we) = electrolyte
    if node.region == Γ_we
        for ia in (nc + 1):(nc + na)
            f[ia] = u[ia]
        end
    end
    return
end


"""
    pnpreaction!(f, u, node, electrolyte)            

Finite volume reaction term
"""
function pnpreaction!(f, u, node, electrolyte)
    (; iϕ, ip) = electrolyte
    f[iϕ] = -chargedensity(u, electrolyte)
    if solvepressure(electrolyte)
        f[ip] = 0
    else
        f[ip] = u[ip]
    end
    return
end

function pnpreaction!(f, u, node, celldata::AbstractCellData)
    elys = electrolytes(celldata)
    pnpreaction!(f, u, node, elys[node.region])
    return
end

"""
     dμex(γk, γl, electrolyte)

Calculate differences of excess chemical potentials from activity coefficients
"""
@inline function dμex(γk, γl, electrolyte)
    (; RT, rlog) = electrolyte
    return (rlog(γk) - rlog(γl)) * RT
end

"""
    μex_flux!(f,dϕ,ck,cl,γk,γl,electrolyte; evelo=0.0)

Excess chemical potential (Sedan) upweind flux, based on modified
Scharfetter-Gummel scheme, see Gaudeul/Fuhrmann 2022.

Apparantly first described by Yu, Zhiping  and Dutton, Robert, SEDAN III, www-tcad.stanford.edu/tcad/programs/sedan3.html

see also the 198? Fortran code available via
https://web.archive.org/web/20210518233152/http://www-tcad.stanford.edu/tcad/programs/oldftpable.html

Verification calculation is in the paper.
"""
function μex_flux!(f, dϕ, ck, cl, γk, γl, electrolyte; evelo = 0.0)
    (; D, z, F, RT, cspecies) = electrolyte
    for ic in cspecies
        bp, bm = fbernoulli_pm(z[ic] * dϕ * F / RT + dμex(γk[ic], γl[ic], electrolyte) / RT - evelo / D[ic])
        f[ic] = D[ic] * (bm * ck[ic] - bp * cl[ic])
    end
    return nothing
end

"""
    act_flux!(ic,dϕ,ck,cl,γk,γl,electrolyte; evelo=0)

Scharfetter-Gummel upwind flux expression based on  activities, see Fuhrmann, CPC 2015.
As shown in Chainais-Hilliaret, Cances, Fuhrmann, Gaudeul, 2019, this is
consistent to thermodynamic equilibrium but shows inferior convergence behavior.
"""
function act_flux!(f, dϕ, ck, cl, γk, γl, electrolyte; evelo = 0.0)
    (; D, z, F, RT, cspecies) = electrolyte
    for ic in cspecies
        Dx = D[ic] * (1 / γk[ic] + 1 / γl[ic]) / 2
        bp, bm = fbernoulli_pm(z[ic] * dϕ * F / RT - evelo / D[ic])
        f[ic] = Dx * (bm * ck[ic] * γk[ic] - bp * cl[ic] * γl[ic])
    end
    return nothing
end

"""
    cent_flux!(ic,dϕ,ck,cl,γk,γl,electrolyte; evelo = 0)

Flux expression based on central differences, see Gaudeul/Fuhrmann 2022.
As shown in the paper  this is
covergent,  consistent to thermodynamic equilibrium but may show 
inferior convergence behavior.
"""
function cent_flux!(f, dϕ, ck, cl, γk, γl, electrolyte; evelo = 0.0)
    (; D, z, F, RT, cspecies, rlog) = electrolyte
    for ic in cspecies
        μk = rlog(ck[ic]) * RT
        μl = rlog(cl[ic]) * RT
        f[ic] = D[ic] * 0.5 * (ck[ic] + cl[ic]) *
            ((μk - μl + dμex(γk[ic], γl[ic], electrolyte) + z[ic] * F * dϕ) / RT - evelo / D[ic])
    end
    return nothing
end

"""
    pnpflux(f, u, edge, electrolyte)

Finite volume flux. It calls either [`μex_flux!`](@ref), [`cent_flux!`](@ref) or [`act_flux!`](@ref).
"""
function pnpflux!(f, u, edge, electrolyte)
    (;
        nc, ip, iϕ,
        ε_0, ε, ε_dec,
        eneutral, pscale, p_bulk,
        upwindflux!, actcoeff!,
        γk_cache, γl_cache,
    ) = electrolyte

    evelo = edgevelocity(electrolyte, edge.index)

    pk, pl = u[ip, 1] * pscale - p_bulk, u[ip, 2] * pscale - p_bulk
    ϕk, ϕl = u[iϕ, 1], u[iϕ, 2]
    ck, cl = view(u, :, 1), view(u, :, 2)
    qk, ql = chargedensity(ck, electrolyte), chargedensity(cl, electrolyte)
    γk, γl = get_tmp(γk_cache, u), get_tmp(γl_cache, u)
    actcoeff!(γk, ck, pk, electrolyte)
    actcoeff!(γl, cl, pl, electrolyte)

    xmid = MVector{3, Float64}(undef)
    for i in 1:size(edge.coord)[1]
        xmid[i] = 0.5 * (edge[i, 1] + edge[i, 2])
    end

    dϕ = ϕk - ϕl
    f[iϕ] = ε_dec(xmid) * ε * ε_0 * dϕ * !eneutral
    if solvepressure(electrolyte)
        f[ip] = u[ip, 1] - u[ip, 2] + (qk + ql) * dϕ / (2 * pscale)
    end

    upwindflux!(f, dϕ, ck, cl, γk, γl, electrolyte; evelo)

    return
end

function pnpflux!(f, u, edge, celldata::AbstractCellData)
    elys = electrolytes(celldata)
    pnpflux!(f, u, edge, elys[edge.region])
    return
end

function sgflux!(y, u, edge, data)
    (; nc, z, iϕ, ip, ε_0, ε, ε_dec, F, RT, D, eneutral) = data
    dϕ = u[iϕ, 1] - u[iϕ, 2]
    xmid = MVector{3, Float64}(undef)
    for i in 1:size(edge.coord)[1]
        xmid[i] = 0.5 * (edge[i, 1] + edge[i, 2])
    end
    y[iϕ] = ε_dec(xmid) * ε * ε_0 * dϕ * !eneutral


    # Dummy equation for pressure
    y[ip] = u[ip, 1] - u[ip, 2]

    for ic in 1:nc
        bp, bm = fbernoulli_pm(z[ic] * dϕ * F / RT)
        y[ic] = D[ic] * (bm * u[ic, 1] - bp * u[ic, 2]) # Scharfetter-Gummel flux
    end
    return nothing
end


"""
    pseudopotentiostat(f0, u0, sys, data)

VoronoiFVM generic operator implementing an IR-compensated pseudopotentiostat.
"""
function pseudopotentiostat(f0, u0, sys, data)
    f0 .= 0.0
    (; iϕ, ϕ_we, i_ref, Γ_we) = data
    f = reshape(f0, sys)
    u = reshape(u0, sys)
    @views potentialbcondition!(f[:, 1], u[:, 1], (; region = Γ_we), data, u[iϕ, i_ref] + ϕ_we)
    return
end

"""
    pnpbstorage_ohmicdrop!(y, u, node, data)

VoronoiFVM boundary storage callback used in the `:ohmicdrop` compensation mode.
"""
function pnpbstorage_ohmicdrop!(y, u, node, data)
    (; Γ_we, iq, icc) = data
    pnpbstorage!(y, u, node, data)
    if node.region == Γ_we
        y[icc] = u[iq]
    end
    return nothing
end


"""
    ohmicdropcompensation(f0, u0, sys, data)

VoronoiFVM generic operator implementing ohmic-drop IR compensation.

The operator couples three global constraints at the electrode boundary `BP1`:

1. **Helmholtz surface charge** (`ia`):
   Species `iq` tracks the double layer  charge.

2. **Capacitive current** (`icc`): via [`pnpbstorage`](@ref) the storage term
   for `icc` is `∂Q/∂t`, giving `j_C = ∂Q/∂t` as the capacitive current

3. **Electrode potential** (`iϕ`): after evaluating the reaction
"""
function ohmicdropcompensation(f0, u0, sys, data)
    f0 .= 0.0
    (; iϕ, iq, icc, ϕ_we, Γ_we, i_ref, Ru, F, ircompfactor, nv, redoxreaction) = data

    f = reshape(f0, sys)
    u = reshape(u0, sys)

    q = zero(eltype(u))


    f[icc, 1] = u[icc, 1]
    j_C = u[icc, 1]

    @views redoxreaction(
        f[:, 1], u[:, 1],
        nothing, data
    )
    ### todo: replace species 1
    j_F = -f[1, 1] * F


    ϕ_DL = ircompfactor * Ru * (j_F + j_C)

    for i in 1:i_ref
        @views q += chargedensity(u[:, i], data) * nv[i]
    end
    if !(data.C_gap ≈ C_large)
        #        q += data.C_gap * (u[iϕ, 1] - (ϕ_DL + ϕ_we))
    end

    f[iq, 1] = u[iq, 1] - q

    @views potentialbcondition!(
        f[:, 1], u[:, 1], (; region = Γ_we), data,
        ϕ_DL + ϕ_we
    )
    return
end


struct PNPSystem <: AbstractElectrochemicalSystem
    vfvmsys::VoronoiFVM.System
end


"""
    PNPSystem(grid;
             celldata=ElectrolyteData(),
             bcondition=(f, u, n, e)-> nothing,
             reaction=(f, u, n, e)-> nothing,
             kwargs...)

Create VoronoiFVM system for generalized Poisson-Nernst-Planck. Input:
- `grid`: discretization grid
- `celldata`: instance of ElectrolyteData or of subtype of AbstractCellData
- `bcondition`: boundary condition
- `reaction` : reactions of the bulk species
- `kwargs`: Keyword arguments of VoronoiFVM.System
"""
function PNPSystem(
        grid::ExtendableGrid;
        celldata::Union{ElectrolyteData, AbstractCellData} = ElectrolyteData(),
        kwargs...
    )
    return PNPSystem(grid, celldata; kwargs...)
end

function PNPSystem(
        grid::ExtendableGrid,
        celldata::AbstractElectrolyteData;
        flux = pnpflux!,
        bcondition = (f, u, n, e) -> nothing,
        reaction = (f, u, n, e) -> nothing,
        kwargs...
    )
    update_derived!(celldata)
    celldata.nv = ones(num_nodes(grid))
    (; ircompensation, iq, icc, Γ_we) = celldata
    ircompensation ∈ (:none, :pseudopotentiostat, :ohmicdrop) ||
        error("ircompensation must be :none, :pseudopotentiostat, or :ohmicdrop, got :$ircompensation")

    # find iref
    coord = grid[Coordinates]
    x_ref = [celldata.x_ref[i] for i in 1:dim_space(grid)]
    dmin = 1.0e30
    imin = 0
    for i in 1:size(coord, 2)
        d = norm(x_ref - coord[:, i])
        if d < dmin
            dmin = d
            imin = i
        end
    end
    celldata.i_ref = imin

    function _pnpreaction!(f, u, node, electrolyte::AbstractElectrolyteData)
        pnpreaction!(f, u, node, electrolyte)
        reaction(f, u, node, electrolyte)
        return nothing
    end

    species = union(celldata.cspecies, [celldata.ip, celldata.iϕ])

    if ircompensation == :none
        sys = VoronoiFVM.System(
            grid;
            data = celldata,
            flux,
            reaction = _pnpreaction!,
            storage = pnpstorage!,
            bstorage = pnpbstorage!,
            bcondition,
            species,
            kwargs...
        )
    elseif ircompensation == :pseudopotentiostat
        sys = VoronoiFVM.System(
            grid;
            data = celldata,
            flux,
            reaction = _pnpreaction!,
            storage = pnpstorage!,
            bstorage = pnpbstorage!,
            bcondition,
            species,
            generic = pseudopotentiostat,
            kwargs...
        )
    elseif ircompensation == :ohmicdrop
        sys = VoronoiFVM.System(
            grid;
            data = celldata,
            flux,
            reaction = _pnpreaction!,
            storage = pnpstorage!,
            bstorage = pnpbstorage_ohmicdrop!,
            bcondition,
            species,
            generic = ohmicdropcompensation,
            kwargs...
        )
    end

    for ia in celldata.sspecies
        enable_boundary_species!(sys, ia, [celldata.Γ_we])
    end
    if ircompensation == :ohmicdrop
        enable_boundary_species!(sys, iq, [Γ_we])
        enable_boundary_species!(sys, icc, [Γ_we])
        celldata.nv .= nodevolumes(sys)
    end
    return PNPSystem(sys)
end

function PNPSystem(
        grid::ExtendableGrid,
        celldata::AbstractCellData;
        bcondition = (f, u, n, e) -> nothing,
        reaction = (f, u, n, e) -> nothing,
        kwargs...
    )
    update_derived!(celldata)
    function _pnpreactionv!(f, u, node, electrolyte::AbstractElectrolyteData)
        pnpreaction!(f, u, node, electrolyte)
        reaction(f, u, node, electrolyte)
        return nothing
    end

    function _pnpreactionv!(f, u, node, celldata::AbstractCellData)
        elys = electrolytes(celldata)
        return _pnpreactionv!(f, u, node, elys[node.region])
    end
    elys = electrolytes(celldata)
    for ely in elys[1:end]
        @assert ely.ip == pressure_index(celldata)
        @assert ely.iϕ == voltage_index(celldata)
    end
    @assert(length(elys) >= num_cellregions(grid))


    sys = VoronoiFVM.System(
        grid;
        data = celldata,
        flux = pnpflux!,
        reaction = _pnpreactionv!,
        storage = pnpstorage!,
        bcondition,
        kwargs...
    )
    enable_species!(sys; species = pressure_index(celldata))
    enable_species!(sys; species = voltage_index(celldata))
    for region in 1:num_cellregions(grid)
        enable_species!(sys; species = elys[region].cspecies, regions = [region])
    end

    return PNPSystem(sys)
end

"""
    electrolytedata(sys)
Extract electrolyte data from system.
"""
electrolytedata(sys::AbstractElectrochemicalSystem) = sys.vfvmsys.physics.data

"""
    celldata(sys)
Extract celldata from system.
"""
celldata(sys::AbstractElectrochemicalSystem) = sys.vfvmsys.physics.data

Base.@deprecate electrolytedata(sys) celldata(sys)


"""
    unknowns(sys)

Return vector of unknowns initialized with bulk data.
"""
function VoronoiFVM.unknowns(esys::AbstractElectrochemicalSystem)
    sys = esys.vfvmsys
    edata = electrolytedata(esys)
    u = unknowns(sys)
    if isa(edata, ElectrolyteData)
        (; iϕ, ip, cspecies, sspecies, c_bulk, Γ_we) = edata
        @views u[iϕ, :] .= 0
        @views u[ip, :] .= 0
        for ic in cspecies
            @views u[ic, :] .= c_bulk[ic]
        end
        for ia in sspecies
            @views u[ia, :] .= 0
        end
    end
    return u
end
