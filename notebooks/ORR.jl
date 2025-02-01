### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
begin
    import Pkg # hide
    Pkg.activate(joinpath(@__DIR__, "..", "docs")) # hide
    using Revise # hide
    using PlutoUI, HypertextLiteral
    using LiquidElectrolytes
    using Printf
    using LessUnitful
    using ExtendableGrids
    using VoronoiFVM
    using GridVisualize
    using StaticArrays
    using Interpolations
    using Test
    if isdefined(Main,:PlutoRunner)
        import CairoMakie	
   	default_plotter!(CairoMakie)
 	CairoMakie.activate!(type="svg")
    end
end

# ╔═╡ 1f5732a6-c15a-4df0-8927-f1e031643d26
md"""
# ORR.jl

[Pluto source](https://raw.githubusercontent.com/j-fu/LiquidElectrolytes.jl/main/notebooks/ORR.jl)

"""

# ╔═╡ b609ba76-e066-4192-a2fc-97b9e659fa12
md"""
Demonstration of oxygen reduction reaction
"""

# ╔═╡ 5482615c-cb5b-4c44-99c4-19416eabae7f
md"""
```math
  4H^+ + O_2 + 4e^- \leftrightharpoons  2H_2O
```
"""

# ╔═╡ 9eebfe1a-2d14-4b9f-a697-71078be8c6d9
pkgdir(LiquidElectrolytes)

# ╔═╡ 84e05551-7d51-4b2c-88f2-b186ad6a244a
md"""
## Setup
"""

# ╔═╡ 16b9af50-11c5-4bdf-b5d8-8d2d9331b5e9
md"""
### Units
"""

# ╔═╡ 515baffb-2e22-401d-aacd-15971dd4365e
begin
    LessUnitful.@phconstants R N_A e
    @unitfactors nm cm μF mol dm s V
    const F = N_A * e
end;

# ╔═╡ bc896b9c-03a9-4da0-be80-f096766cb039
md"""
### Data
"""

# ╔═╡ 50ff0114-f9b1-4bd0-8190-153fef07ef3a
begin
    const vmin = -1V
    const vmax = 1V
    const vdelta = 0.025 * V
    const molarity = 0.1
    const nref = 0
    const κ = 10.0
    const vfac = 1
    const scheme = :μex
    const R0 = 5.0e-16mol / (cm^2 * s)

    const Δg = 0.0
    const β = 0.5
    const ϕ_we = 0.0
    const ihplus = 1
    const iso4 = 2
    const io2 = 3
    const z = [1, -2, 0]
    const allκ = [κ, κ, 0]
end;

# ╔═╡ 231ebfc2-51b6-4027-8be2-75891dfed7e0
md"""
### Solver control
"""

# ╔═╡ a4e0379d-f7f6-4b61-bf38-5eb17f67505a
solver_control = (max_round = 4,
                  tol_round = 1.0e-8,
                  reltol = 1.0e-8,
                  abstol = 1.0e-9,
                  verbose = "",
                  maxiters = 20)

# ╔═╡ 970389b5-d2c1-4992-9978-aca1ccd3d2fc
md"""
### Reaction description
"""

# ╔═╡ 136cf1aa-75c6-4aaa-b93b-801391ec800c
md"""
In the following reaction function, the  balance with the solvent is fulfilled automatically via the incompressibility constraint. Any material which is removed via the boundary reaction is automatically replaced by solvent with the corresponding volume. So "ignoring the solvent" here is correct.
"""

# ╔═╡ b916eb92-8ba8-49aa-bd7c-1bfc91e813d4
function halfcellbc(f, u, bnode, data)
    bulkbcondition(f, u, bnode, data)
    (; iϕ, eneutral, ϕ_we, Γ_we, RT) = data

    if bnode.region == Γ_we
        f .= 0.0
        if !eneutral
            boundary_dirichlet!(f, u, bnode; species = iϕ, region = Γ_we,
                                value = ϕ_we)
        end
        μh2o, μ = chemical_potentials!(MVector{4, eltype(u)}(undef), u, data)
        A = (4 * μ[ihplus] + μ[io2] - 2μh2o + Δg +
             4 * eneutral * F * (u[iϕ] - ϕ_we)) / (RT)
        r = rrate(R0, β, A)
        f[ihplus] -= 4 * r
        f[io2] -= r
    end
end

# ╔═╡ 420a82e0-5fc2-47ea-8916-d88910655d50
md"""
## Nernst-Planck halfcell
"""

# ╔═╡ 392a648c-12a2-41fb-b3f6-aa3cfe3cbcd7
grid = let
    hmin = 1.0e-1 * nm * 2.0^(-nref)
    hmax = 1.0 * nm * 2.0^(-nref)
    L = 20.0 * nm
    X = geomspace(0, L, hmin, hmax)
    simplexgrid(X)
end

# ╔═╡ 60d85155-9fa7-4740-9769-212ceef1918b
begin
    celldata = ElectrolyteData(;
                               nc = 3,
                               z,
                               κ = allκ,
                               Γ_we = 1,
                               Γ_bulk = 2,
                               eneutral = false,
	                            log=RLog(),
                               scheme)

    celldata.v .*= vfac
    celldata.weights[celldata.ip] = 0
    celldata.weights[1:(celldata.nc)] .= 0
    c_bulk = celldata.c_bulk
    c_bulk[io2] = 0.001 * mol / dm^3
    c_bulk[iso4] = molarity * mol / dm^3
    c_bulk[ihplus] = 2.0 * molarity * mol / dm^3
end

# ╔═╡ 36306b3d-681f-423c-9b32-db6562c5c157
@test isincompressible(celldata.c_bulk, celldata)

# ╔═╡ 6e8f10b8-617e-4879-a6cb-fe93a4fd7226
@test iselectroneutral(celldata.c_bulk, celldata)

# ╔═╡ 1d23b6b5-88cf-4200-b612-82dffcc5cca7
cell = PNPSystem(grid; bcondition = halfcellbc, celldata)

# ╔═╡ 763c393c-e0c8-447e-92e4-f2a5f0de2a30
result=LiquidElectrolytes.ivsweep(cell;store_solutions=true,
                                                voltages = vmin:vdelta:vmax,
                                                solver_control...)

# ╔═╡ 22bc5f42-1a21-41a5-a059-b4ac44a29566
md"""
### Results
"""

# ╔═╡ d7b10140-7db7-4be0-88c3-53ba1f203310
@bind vshow PlutoUI.Slider(range(extrema(result.voltages)...; length = 101),
                           show_value = true)

# ╔═╡ 300ed474-76c5-47e9-b15a-8c4c93082268
md"""
### Plotting functions
"""

# ╔═╡ f1857d7d-cec5-42a5-88d6-1d1f620f894c
curr(J, ix) = [F * j[ix] for j in J]

# ╔═╡ 5bc4f11f-24c6-4af8-a554-1b5771f1f2b0
function curr_h2o(J)
    -4κ * curr(J, ihplus) - (κ + 2) * curr(J, io2)
end

# ╔═╡ b81676e8-dcec-49fd-b350-f26ac61243ec
function plotcurr(result)
    scale = 1 / (mol / dm^3)
    volts = result.voltages
    vis = GridVisualizer(;
                         size = (600, 300),
                         tilte = "IV Curve",
                         xlabel = "Φ_WE/V",
                         ylabel = "I",
                         legend = :lt)
    scalarplot!(vis,
                volts,
                curr(result.j_bulk, ihplus);
                linestyle = :dash,
                label = "H+, bulk",
                color = :red)
    scalarplot!(vis,
                volts,
                curr(result.j_we, ihplus);
                color = :red,
                clear = false,
                linestyle = :solid,
                label = "H+, we")

    scalarplot!(vis,
                volts,
                curr(result.j_bulk, io2);
                linestyle = :dash,
                label = "O2, bulk",
                color = :green,
                clear = false)
    scalarplot!(vis,
                volts,
                curr(result.j_we, io2);
                color = :green,
                clear = false,
                linestyle = :solid,
                label = "O2, we")

    scalarplot!(vis,
                volts,
                curr(result.j_bulk, io2);
                linestyle = :dash,
                label = "O2, bulk",
                color = :green,
                clear = false)
    scalarplot!(vis,
                volts,
                curr(result.j_we, io2);
                color = :green,
                clear = false,
                linestyle = :solid,
                label = "O2, we")

    scalarplot!(vis,
                volts,
                curr_h2o(result.j_bulk) / 100;
                linestyle = :dash,
                label = "H2O/100, bulk",
                color = :blue,
                clear = false)
    scalarplot!(vis,
                volts,
                curr_h2o(result.j_we) / 100;
                color = :blue,
                clear = false,
                linestyle = :solid,
                label = "H2O/100, we")

    reveal(vis)
end

# ╔═╡ 7891a252-8fdf-40df-a205-64ca4078a542
plotcurr(result)

# ╔═╡ e5ea40e3-ef3f-4c89-bc69-91e95dfb3cd0
function plot1d!(vis,result, celldata, vshow)
    vinter = linear_interpolation(result.voltages, [j[io2] for j in result.j_we])
    tsol=LiquidElectrolytes.voltages_solutions(result)

               sol = tsol(vshow)
            c0 = solventconcentration(sol, celldata)
            scale = 1.0 / (mol / dm^3)
	    ishow=vinter(vshow)
            title = "Φ_we=$(round(vshow,digits=4)), I=$(round(vinter(vshow),sigdigits=4))"
	#    title = @sprintf("Φ_we=%+1.2f I=%+1.4f",vshow,ishow)
            
            scalarplot!(vis, grid, sol[io2, :] * scale.|>abs; color = :green, label = "O_2",title)
            scalarplot!(vis,
                        grid,
                        sol[iso4, :] * scale .|>abs;
                        color = :gray,
                        clear = false,
                        label = "SO4--")
            scalarplot!(vis,
                        grid,
                        sol[ihplus, :] * scale .|>abs;
                        color = :red,
                        clear = false,
                        label = "H+")
            scalarplot!(vis, grid, c0 * scale.|>abs; color = :blue, clear = false,
                        label = "H2O")	
end

# ╔═╡ 9226027b-725d-446e-bc14-dd335a60ec09
function plot1d(result,celldata, vshow)
    vis = GridVisualizer(;
                         size = (600, 250),
                         yscale = :log,
                         limits = (1.0e-6, 100),
                         legend = :rt)
    
	   
     plot1d!(vis,result,celldata,vshow)
     reveal(vis)
end

# ╔═╡ cc80544f-ca62-49fb-b907-4bd194b11ee5
function plot1d(result,celldata)
    vinter = linear_interpolation(result.voltages, [j[io2] for j in result.j_we])
    tsol=LiquidElectrolytes.voltages_solutions(result)
    vis = GridVisualizer(;
                         size = (600, 250),
                         yscale = :log,
                         limits = (1.0e-6, 100),
                         legend = :rt)
    
    video="orr.gif"
    vrange=range(extrema(result.voltages)...; length = 101)
    
    movie(vis,file="orr.mp4") do vis
   	for vshow in vrange
    plot1d!(vis,result,celldata,vshow)
		reveal(vis)
	end
    end
	isdefined(Main,:PlutoRunner) && LocalResource("orr.mp4")
end

# ╔═╡ 56eb52b1-9017-4485-83d6-b7ef15ad522f
plot1d(result, celldata, vshow)

# ╔═╡ 3f639537-21e9-4dc0-8eeb-59e7d28afee1
plot1d(result, celldata)

# ╔═╡ 1ac7646a-76ae-4e8f-9d9d-ecaccc262857
function cplot(cell, result)
    scale = 1.0 / (mol / dm^3)
    tsol=LiquidElectrolytes.voltages_solutions(result)
    j_we=result.j_we
    currs = curr(j_we, io2)
    vis = GridVisualizer(; resolution = (1200, 400), layout = (1, 3),
                         gridscale = 1.0e9)
    xmax = 10 * nm
    xlimits = [0, xmax]
    aspect = 2 * xmax / (tsol.t[end] - tsol.t[begin])
    scalarplot!(vis[1, 1],
                cell,
                tsol;
                scale,
                species = io2,
                aspect,
                xlimits,
                title = "O2",
                colormap = :summer)
    scalarplot!(vis[1, 2],
                cell,
                tsol;
                species = ihplus,
                aspect,
                scale,
                xlimits,
                title = "H+",
                colormap = :summer)
    scalarplot!(vis[1, 3],
                1000 * tsol[io2, 1, :] * scale,
                tsol.t;
                label = "1000*O2",
                xlabel = "c",
                color = :green,
                clear = false)
    scalarplot!(vis[1, 3],
                tsol[ihplus, 1, :] * scale,
                tsol.t;
                title = "c(0)",
                xlabel = "c",
                ylabel = "V",
                label = "H+",
                color = :red,
                clear = false)
    scalarplot!(vis[1, 3],
                tsol[iso4, 1, :] * scale,
                tsol.t;
                label = "SO4--",
                color = :blue,
                clear = false,
                legend = :ct)
    reveal(vis)
end

# ╔═╡ 556c47ee-e172-483b-b922-a6422a0c405f
cplot(cell, result)

# ╔═╡ 1317a982-c416-4d44-804a-8694cc2bbef2
md"""
## Electroneutral case
"""

# ╔═╡ 950f43ba-6555-463a-bed7-36511e17e882
begin
    ncelldata = deepcopy(celldata)
    ncelldata.eneutral = true
    ncell = PNPSystem(grid; bcondition = halfcellbc, celldata = ncelldata)
end

# ╔═╡ 33ee0ded-5bc8-4fe7-bd2f-1cc44bc73f78

nresult = LiquidElectrolytes.ivsweep(ncell;
                                                store_solutions=true,
                                                   voltages = vmin:vdelta:vmax,
                                                   solver_control...)

# ╔═╡ 6a0ff3ea-25af-4682-a8f6-40c481b53d8d
md"""
### Results
"""

# ╔═╡ 63dd0cef-7acd-4507-bbc6-3976181a143d
plotcurr(nresult)

# ╔═╡ c7185947-56ea-4e79-a619-03bf77d5219d
@bind nvshow PlutoUI.Slider(range(extrema(nresult.voltages)...; length = 101),
                            show_value = true)

# ╔═╡ c6b9f3ce-dd3e-474e-b947-3daacc5cd1d0
plot1d(nresult,ncelldata, nvshow)

# ╔═╡ b729b190-a7ed-48b9-9584-fe5271e5dfa4
cplot(ncell, nresult)

# ╔═╡ 42dda2f6-ea60-4cbc-8372-fafd4a1218a8
md"""
### Comparison
"""

# ╔═╡ f4ad49b2-be13-4f6b-97ce-d24eae913279
let vis = GridVisualizer(;
                         size = (600, 400),
                         title = "IV Curve",
                         xlabel = "Φ_WE/V",
                         ylabel = "I/(A/m^2)",
                         legend = :lt,
                         )
    scalarplot!(vis,
                nresult.voltages,
                currents(nresult,io2);
                label = "O2,electroneutral",
                color = :green,
                clear = false,
				limits= [-5.0e4, 5.0e4])

    scalarplot!(vis, result.voltages, currents(result, io2); label = "O2,PNP", color = :red,
                clear = false)
    reveal(vis)
end

# ╔═╡ f9b4d4dc-7def-409f-b40a-f4eba1163741
TableOfContents()

# ╔═╡ 7a93e9a8-8a2d-4b11-84ef-691706c0eb0f
begin
    hrule() = html"""<hr>"""
    function highlight(mdstring, color)
        htl"""<blockquote style="padding: 10px; background-color: $(color);">$(mdstring)</blockquote>"""
    end

    macro important_str(s)
        :(highlight(Markdown.parse($s), "#ffcccc"))
    end
    macro definition_str(s)
        :(highlight(Markdown.parse($s), "#ccccff"))
    end
    macro statement_str(s)
        :(highlight(Markdown.parse($s), "#ccffcc"))
    end

    html"""
        <style>
         h1{background-color:#dddddd;  padding: 10px;}
         h2{background-color:#e7e7e7;  padding: 10px;}
         h3{background-color:#eeeeee;  padding: 10px;}
         h4{background-color:#f7f7f7;  padding: 10px;}
        
	     pluto-log-dot-sizer  { max-width: 655px;}
         pluto-log-dot.Stdout { background: #002000;
	                            color: #10f080;
                                border: 6px solid #b7b7b7;
                                min-width: 18em;
                                max-height: 300px;
                                width: 675px;
                                    overflow: auto;
 	                           }
	
    </style>
"""
end

# ╔═╡ 5beb3a0d-e57a-4aea-b7a0-59b8ce9ff5ce
hrule()


# ╔═╡ Cell order:
# ╟─1f5732a6-c15a-4df0-8927-f1e031643d26
# ╟─b609ba76-e066-4192-a2fc-97b9e659fa12
# ╟─5482615c-cb5b-4c44-99c4-19416eabae7f
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╠═9eebfe1a-2d14-4b9f-a697-71078be8c6d9
# ╟─84e05551-7d51-4b2c-88f2-b186ad6a244a
# ╟─16b9af50-11c5-4bdf-b5d8-8d2d9331b5e9
# ╠═515baffb-2e22-401d-aacd-15971dd4365e
# ╟─bc896b9c-03a9-4da0-be80-f096766cb039
# ╠═50ff0114-f9b1-4bd0-8190-153fef07ef3a
# ╟─231ebfc2-51b6-4027-8be2-75891dfed7e0
# ╠═a4e0379d-f7f6-4b61-bf38-5eb17f67505a
# ╟─970389b5-d2c1-4992-9978-aca1ccd3d2fc
# ╟─136cf1aa-75c6-4aaa-b93b-801391ec800c
# ╠═b916eb92-8ba8-49aa-bd7c-1bfc91e813d4
# ╟─420a82e0-5fc2-47ea-8916-d88910655d50
# ╠═392a648c-12a2-41fb-b3f6-aa3cfe3cbcd7
# ╠═60d85155-9fa7-4740-9769-212ceef1918b
# ╠═36306b3d-681f-423c-9b32-db6562c5c157
# ╠═6e8f10b8-617e-4879-a6cb-fe93a4fd7226
# ╠═1d23b6b5-88cf-4200-b612-82dffcc5cca7
# ╠═763c393c-e0c8-447e-92e4-f2a5f0de2a30
# ╟─22bc5f42-1a21-41a5-a059-b4ac44a29566
# ╟─7891a252-8fdf-40df-a205-64ca4078a542
# ╠═d7b10140-7db7-4be0-88c3-53ba1f203310
# ╠═56eb52b1-9017-4485-83d6-b7ef15ad522f
# ╠═3f639537-21e9-4dc0-8eeb-59e7d28afee1
# ╟─556c47ee-e172-483b-b922-a6422a0c405f
# ╟─300ed474-76c5-47e9-b15a-8c4c93082268
# ╟─f1857d7d-cec5-42a5-88d6-1d1f620f894c
# ╟─5bc4f11f-24c6-4af8-a554-1b5771f1f2b0
# ╟─b81676e8-dcec-49fd-b350-f26ac61243ec
# ╠═e5ea40e3-ef3f-4c89-bc69-91e95dfb3cd0
# ╠═9226027b-725d-446e-bc14-dd335a60ec09
# ╠═cc80544f-ca62-49fb-b907-4bd194b11ee5
# ╟─1ac7646a-76ae-4e8f-9d9d-ecaccc262857
# ╟─1317a982-c416-4d44-804a-8694cc2bbef2
# ╠═950f43ba-6555-463a-bed7-36511e17e882
# ╠═33ee0ded-5bc8-4fe7-bd2f-1cc44bc73f78
# ╟─6a0ff3ea-25af-4682-a8f6-40c481b53d8d
# ╠═63dd0cef-7acd-4507-bbc6-3976181a143d
# ╠═c7185947-56ea-4e79-a619-03bf77d5219d
# ╟─c6b9f3ce-dd3e-474e-b947-3daacc5cd1d0
# ╟─b729b190-a7ed-48b9-9584-fe5271e5dfa4
# ╟─42dda2f6-ea60-4cbc-8372-fafd4a1218a8
# ╠═f4ad49b2-be13-4f6b-97ce-d24eae913279
# ╟─5beb3a0d-e57a-4aea-b7a0-59b8ce9ff5ce
# ╟─f9b4d4dc-7def-409f-b40a-f4eba1163741
# ╟─7a93e9a8-8a2d-4b11-84ef-691706c0eb0f
