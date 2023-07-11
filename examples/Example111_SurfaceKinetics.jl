#=
# Suface kinetics
([source code](SOURCE_URL))

=#

module Example111_SurfaceKinetics
using LessUnitful
using ExtendableGrids, GridVisualize
using VoronoiFVM
using LiquidElectrolytes
using Colors
using StaticArrays
using InteractiveUtils
using ForwardDiff


function main(;
    nref = 0,
    voltages = (-3.0:0.1:3.0) * ufac"V",
    molarities = [0.001, 0.01, 0.1, 1],
    scheme = :μex,
    κ = 10.0,
    Plotter = nothing,
    new = false,
    kwargs...,
)

    @local_phconstants N_A e R ε_0 k_B
    F = N_A * e
    c_0 = 2.99792458e8
    @local_unitfactors cm μF mol dm s mA A nm bar eV μA



    defaults = (;
        max_round = 3,
        tol_round = 1.0e-9,
        verbose = "e",
        reltol = 1.0e-8,
        tol_mono = 1.0e-10,
    )

    kwargs = merge(defaults, kwargs)

    hmin = 1.0e-1 * ufac"μm" * 2.0^(-nref)
    hmax = 1.0 * ufac"μm" * 2.0^(-nref)
    L = 80.0 * ufac"μm"
    X = geomspace(0, L, hmin, hmax)
    grid = simplexgrid(X)

    T = 273.15 + 25 * ufac"K"


    ## kinetic model
    ## A+_aq <-> A+_ads,	                    ##1
    ## A+_ads + e- <-> A_ads,                ##2
    ## A_ads <-> A_aq                        ##3


    bulk_species = ["A+_aq", "A_aq"]
    surface_species = ["A+_ads", "A_ads"]

    C_gap = 20 * ufac"μF/cm^2"
    ϕ_pzc = 0.2 * ufac"V"

    iaplus = 1
    ibminus = 2
    ia = 3
    iaplus_ads = 4
    ia_ads = 5


    function halfcellbc(
        f,
        u::VoronoiFVM.BNodeUnknowns{Tval,Tv,Tc,Tp,Ti},
        bnode,
        data,
    ) where {Tval,Tv,Tc,Tp,Ti}
        (; nc, na, Γ_we, Γ_bulk, ϕ_we, ip, iϕ, v, v0, T, RT, ε) = data

        bulkbcondition(f, u, bnode, data; region = Γ_bulk)


        if bnode.region == Γ_we

            ## boundary_dirichlet!(f,u,bnode;species=iϕ,region=Γ_we,value=ϕ_we)

            ## Robin b.c. for the Poisson equation

            boundary_robin!(f, u, bnode, iϕ, C_gap / ε, C_gap * (ϕ_we - ϕ_pzc) / ε)

            ## surface current density
            sigma = C_gap * (ϕ_we - u[iϕ] - ϕ_pzc)

            kf = zeros(Tval, 3)
            kr = zeros(Tval, 3)

            ## A+_aq <-> A+_ads,	                    ##1
            ΔG_ads_aplus = 1.0 * eV + 1.0e+1 * sigma * eV

            kf[1] = 1.0e13 * exp(-max(ΔG_ads_aplus, 0.0) / (k_B * T))
            kr[1] = 1.0e13 * exp(-max(-ΔG_ads_aplus, 0.0) / (k_B * T))

            ## A+_ads + e- <-> A_ads,                ##2          
            ΔG_rxn = 1.0 * eV + 1.5e+1 * sigma * eV

            kf[2] = 1.0e13 * exp(-max(ΔG_rxn, 0.0) / (k_B * T))
            kr[2] = 1.0e13 * exp(-max(-ΔG_rxn, 0.0) / (k_B * T))

            ## A_ads <-> A_aq                        ##3
            ΔG_ads_a = 0.5 * eV + 0.5e+1 * sigma * eV

            kf[3] = 1.0e13 * exp(-max(ΔG_ads_a, 0.0) / (k_B * T))
            kr[3] = 1.0e13 * exp(-max(-ΔG_ads_a, 0.0) / (k_B * T))

            S = 1.0e-6 / N_A * (1.0e10)^2 * ufac"mol/m^2"
            θ_free = 1 - u[iaplus_ads] - u[ia_ads]

            rates = zeros(Tval, 3)

            rates[1] = kf[1] * u[iaplus] * θ_free - kr[1] * u[iaplus_ads]
            rates[2] = kf[2] * u[iaplus_ads] - kr[2] * u[ia_ads]
            rates[3] = kf[3] * u[ia_ads] - kr[3] * u[ia] * θ_free

            println(
                "rate constants: $(ForwardDiff.value.(kf)) and $(ForwardDiff.value.(kr))",
            )
            println("rates: $(ForwardDiff.value.(rates))")
            println("aplus: $(ForwardDiff.value(u[iaplus]))")

            ## bulk species
            f[iaplus] += -rates[1] * S
            f[ia] += rates[3] * S

            ## surface species
            f[iaplus_ads] += rates[1] - rates[2]
            f[ia_ads] += rates[2] - rates[3]

        end
        nothing
    end


    celldata = ElectrolyteData(;
        nc = 3,
        na = 2,
        z = [1, -1, 0],
        D = [2.0e-9, 2.0e-9, 2.0e-9] * ufac"m^2/s", ## from Ringe paper
        T = T,
        eneutral = false,
        κ = fill(κ, 3),
        Γ_we = 1,
        Γ_bulk = 2,
        scheme,
    )

    (; iϕ::Int, ip::Int) = celldata

    celldata.c_bulk[iaplus] = 0.1 * mol / dm^3
    celldata.c_bulk[ibminus] = 0.1 * mol / dm^3
    celldata.c_bulk[ia] = 0.1 * mol / dm^3

    @assert isapprox(celldata.c_bulk' * celldata.z, 0, atol = 1.0e-10)

    cell = PNPSystem(grid; bcondition = halfcellbc, celldata)

    un = unknowns(cell)
    @views un[iaplus_ads, :] .= 0.1
    @views un[ia_ads, :] .= 0.1


    result = ivsweep(cell; voltages, store_solutions=true, kwargs...)
    tsol = voltages_solutions(result)
    currs = currents(result, ia)
    volts = result.voltages
    
    vis = GridVisualizer(Plotter = Plotter, layout = (1, 1))

    scalarplot!(
        vis[1, 1],
        volts,
        currs * ufac"cm^2/mA",
        color = :red,
        markershape = :utriangle,
        markersize = 7,
        markevery = 10,
        label = "PNP",
        legend = :lt,
        xlabel = "Δϕ/V",
        ylabel = "I/(mA/cm^2)",
    )
    ##scalarplot!(vis[2,1], sigmas, energies, color="black",clear=true,xlabel="σ/(μC/cm^s)",ylabel="ΔE/eV")
    ##scalarplot!(vis[2,1], ϕs, rs, xlimits=(-1.5,-0.6), yscale=:log, xlabel="Δϕ/V", ylabel="c(CO2)/M")
    for curr in currs
        println(curr)
    end
    return reveal(vis)
end

end
