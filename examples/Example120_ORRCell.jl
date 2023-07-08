#=
# ORR Half cell
([source code](SOURCE_URL))

I-V sweep for Oxygen Reduction

Methods called:
- [`ElectrolyteData`](@@ref)
- [`ivsweep`](@@ref)
- [`dlcapsweep`](@@ref)
- [`PNPSystem`](@@ref)
=#


module Example120_ORRCell
using ExtendableGrids, GridVisualize
using VoronoiFVM
using LiquidElectrolytes
using Colors
using StaticArrays
using LessUnitful


function main(;
              voltages = -1:0.1:1,
              compare = false,
              molarity = 0.1,
              nref = 0,
              κ = 10.0,
              vfac = 1.0,
              eneutral = false,
              scheme = :μex,
              Plotter = nothing,
              R0::Float64 = 4.0e-15,
              epsreg = 1.0e-20,
              kwargs...,
              )

    @local_phconstants R N_A e
    @local_unitfactors nm cm μF mol dm s
    F = N_A * e


    defaults = (;
        max_round = 3,
        tol_round = 1.0e-10,
        reltol = 1.0e-7,
        tol_mono = 1.0e-7,
        verbose = "e",
    )
    kwargs = merge(defaults, kwargs)

    hmin = 1.0e-1 * nm * 2.0^(-nref)
    hmax = 1.0 * nm * 2.0^(-nref)
    L = 20.0 * nm
    X = geomspace(0, L, hmin, hmax)
    grid = simplexgrid(X)


    R0 = R0 * mol / (cm^2 * s)
    Δg = 0.0
    β = 0.5
    ϕ_we = 0.0
    ihplus = 1
    iso4 = 2
    io2 = 3
    z = [1, -2, 0]
    κ = [κ, κ, 0]

    
    function halfcellbc(f, u, bnode, data)
        bulkbcondition(f, u, bnode, data)
        (; iϕ, eneutral, ϕ_we, Γ_we, RT) = data

        if bnode.region == Γ_we
            f .= 0.0
            if !data.eneutral
                boundary_dirichlet!(
                    f,
                    u,
                    bnode;
                    species = iϕ,
                    region = data.Γ_we,
                    value = data.ϕ_we,
                )
            end
            μh2o, μ = chemical_potentials!(MVector{4,eltype(u)}(undef), u, data)
            A =
                (4 * μ[ihplus] + μ[io2] - 2μh2o + Δg + 4*eneutral * F * (u[iϕ] - data.ϕ_we)) /
                (RT)

            r = rrate(R0, β, A)
            f[ihplus] -= 4 * r
            f[io2] -= r
        end
    end




    celldata =
        ElectrolyteData(; nc = 3, z, κ, Γ_we = 1, Γ_bulk = 2, eneutral, scheme, epsreg)

    celldata.v*=vfac
    
    (; iϕ, c_bulk) = celldata


    c_bulk[io2] = 0.001 * mol / dm^3
    c_bulk[iso4] = molarity * mol / dm^3
    c_bulk[ihplus] = 2.0 * molarity * mol / dm^3


    @assert isapprox(celldata.c_bulk' * celldata.z, 0, atol = 1.0e-12)

    cell = PNPSystem(grid; bcondition = halfcellbc, celldata)


    ## Compare electroneutral and double layer cases
    if compare
        celldata.eneutral = false
        sol=ivsweep(cell; voltages, kwargs...)
        currs=currents(sol,io2)
        volts=sol.voltages
        
        celldata.eneutral = true
        nsol  = ivsweep(cell; voltages, kwargs...)
        ncurrs=currents(nsol,io2)
        nvolts=nsol.voltages
        
        vis = GridVisualizer(;
            Plotter,
            resolution = (600, 400),
            clear = true,
            legend = :lt,
            xlabel = "Δϕ/V",
            ylabel = "I/(A/m^2)",
        )
        scalarplot!(
            vis,
            volts,
            currs,
            color = "red",
            markershape = :utriangle,
            markersize = 7,
            markevery = 10,
            label = "PNP",
        )
        scalarplot!(
            vis,
            nvolts,
            ncurrs,
            clear = false,
            color = :green,
            markershape = :none,
            label = "NNP",
        )
        return reveal(vis)
    end

    ## IVsweep
    vis = GridVisualizer(;Plotter, resolution = (1000, 300), layout = (1, 4))
    
    result = ivsweep(cell; voltages,store_solutions=true, kwargs...)
    currs=LiquidElectrolytes.currents(result,io2)
    sol=LiquidElectrolytes.voltage_solution(result)
    volts=result.voltages

    xmax = 10 * nm
    xlimits = [0, xmax]
    aspect = 2 * xmax / (volts[end] - volts[begin])

    
    scalarplot!(
        vis[1, 1],
        currs,
        volts,
        markershape = :none,
        title = "IV",
        xlabel = "I",
        ylabel = "V",
    )
    scalarplot!(
        vis[1, 2],
        cell,
        sol;
        scale= 1.0/(mol/dm^3),
        species = io2,
        aspect,
        xlimits,
        title = "O2",
        colormap = :summer,
    )
    scalarplot!(
        vis[1, 3],
        cell,
        sol;
        species = ihplus,
        aspect,
        scale= 1.0/(mol/dm^3),
        xlimits,
        title = "H+",
        colormap = :summer,
    )

    scalarplot!(
        vis[1, 4],
        sol[io2, 1, :]*1000,
        volts,
        xlabel = "c",
        label = "1000 O2",
        color = :green,
        clear = false,
        legend = :rc,
    )
    scalarplot!(
        vis[1, 4],
        sol[ihplus, 1, :],
        volts,
        title = "c(0)",
        xlabel = "c",
        ylabel = "V",
        label = "H+",
        color = :red,
        clear=false,
    )
    scalarplot!(
        vis[1, 4],
        sol[iso4, 1, :],
        volts,
        label = "SO4--",
        color = :blue,
        clear = false,
        legend = :rc,
    )

    reveal(vis)
end

end


#=
```@example Example120_ORRCell_1
using Example120_ORRCell,CairoMakie #hide
CairoMakie.activate!(type="svg",visible=false) # hide
Example120_ORRCell.main(Plotter=CairoMakie)
```
=#



#=
```@example Example120_ORRCell_2
using Example120_ORRCell,CairoMakie #hide
CairoMakie.activate!(type="svg",visible=false) # hide
Example120_ORRCell.main(compare=true,Plotter=CairoMakie)
```
=#

