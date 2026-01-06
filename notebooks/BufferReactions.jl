### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ dbaff000-beeb-46c3-bbab-7c7e3a3e9791
# ╠═╡ skip_as_script = true
#=╠═╡
begin
    import Pkg  # hide
    Pkg.activate(joinpath(@__DIR__, "..", "docs")) # hide
    using Revise # hide
	using GridVisualize
	     using CairoMakie	
	   	default_plotter!(CairoMakie)
	 	CairoMakie.activate!(type="svg")
   
end
  ╠═╡ =#

# ╔═╡ 3b5a64e6-24cd-423e-aad8-8f400b338867
begin
    using LessUnitful
    using ExtendableGrids
    using LiquidElectrolytes
    using VoronoiFVM
    using InteractiveUtils
    using ForwardDiff
    using PlutoUI, HypertextLiteral
    using DataFrames
    using Test
    using Printf
end

# ╔═╡ eea5aa29-85d4-4644-bfe4-4b5244eeb045
md"""
# BufferReactions.jl

[Pluto-source](https://raw.githubusercontent.com/j-fu/LiquidElectrolytes.jl/main/notebooks/BufferReactions.jl)
"""

# ╔═╡ 8d717a52-d9d5-491b-8ce8-98da978bbca5
md"""
Demonstration of buffer reactions within the electrolyte of a half-cell in thermodynamical equilibrium
"""

# ╔═╡ b263c6bf-3606-4c3d-9fe1-5d074b294519
pkgdir(LiquidElectrolytes)

# ╔═╡ 4c355a5e-ccb8-4d23-be5c-adc1aa5f8c96
md"""
## Setup
"""

# ╔═╡ 327ed52b-6c81-4c2f-8b79-824017351020
md"""
### Units
"""

# ╔═╡ e8de40fe-472e-437f-92c8-5f75f5d58601
begin
    @unitfactors mol dm eV μA μF cm μm K V m s
    @phconstants c_0 N_A e k_B
    F = N_A * e
end;

# ╔═╡ 28309673-e34e-4bf3-9c53-d04a283364ef
md"""
### Data
"""

# ╔═╡ 2a787207-55ab-4e5b-a9f2-22261ee28a42
begin
    const pH = 6.8
    const T = (273.15 + 25) * K
    const C_gap = 20 * μF / cm^2
    const ϕ_pzc = 0.16 * V
    const vmin = -1.5 * V
    const vmax = 0.0 * V
    const vdelta = 0.1 * V
    const nref = 0

    const ikplus = 1
    const ihco3 = 2
    const ico3 = 3
    const ico2 = 4
    const ico = 5
    const iohminus = 6
    const ihplus = 7
    const nc = 7

    ## CO2 + OH- <=> HCO3-
    const kbe1 = 4.44e7 / (mol / dm^3)
    const kbf1 = 5.93e3 / (mol / dm^3) / s
    const kbr1 = kbf1 / kbe1
    ## HCO3- + OH- <=> CO3-- + H2O
    const kbe2 = 4.66e3 / (mol / dm^3)
    const kbf2 = 1.0e8 / (mol / dm^3) / s
    const kbr2 = kbf2 / kbe2

    ## CO2 + H20 <=> HCO3- + H+
    const kae1 = 4.44e-7 * (mol / dm^3)
    const kaf1 = 3.7e-2 / s
    const kar1 = kaf1 / kae1
    ## HCO3- <=> CO3-- + H+
    const kae2 = 4.66e-5 / (mol / dm^3)
    const kaf2 = 59.44e3 / (mol / dm^3) / s
    const kar2 = kaf2 / kae2

    ## autoprotolyse
    const kwe = 1.0e-14 * (mol / dm^3)^2
    const kwf = 2.4e-5 * (mol / dm^3) / s
    const kwr = kwf / kwe
end;

# ╔═╡ 960fe09f-8bfe-4773-a57e-cbd7774050ce
const bulk = DataFrame(
    :name => ["K⁺", "HCO₃⁻", "CO₃²⁻", "CO₂", "CO", "OH⁻", "H⁺"],
    :z => [1, -1, -2, 0, 0, -1, 1],
    :D => [1.957e-9, 1.185e-9, 0.923e-9, 1.91e-9, 2.23e-9, 5.273e-9, 9.31e-9] * m^2 / s,
    :κ => [8.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0],
    :c_bulk => [0.09105350460641519, 0.091, 2.68e-5, 0.033, 0.0, 10^(pH - 14), 10^(-pH)] * mol / dm^3
)

# ╔═╡ e969bc7f-f0f3-4f70-a522-133ee56cce5e
md"""
### Buffer reactions

Mean-field approach using the law of mass action
"""

# ╔═╡ 3161211b-3188-4e4f-a12e-5a16cfd5e932
md"""
Bicarbonate buffer system in base:

$CO_2 + OH^- \rightleftharpoons HCO_3^-$ 

$HCO_3^- + OH^- \rightleftharpoons CO_3^{2-} + H_2O$
"""

# ╔═╡ 18e13140-2693-4a2b-b8e2-24ae437414cd
md"""
Bicarbonate buffer system in acid

$CO_2 + H_2O \rightleftharpoons HCO_3^- + H^+$ 

$HCO_3^- \rightleftharpoons CO_3^{2-} + H^+$
"""

# ╔═╡ c1892a61-dd26-4173-8675-8136c8360469
md"""
Autoprotolysis of Water

$H_2O \rightleftharpoons H^+ + OH^-$
"""

# ╔═╡ 12c99cb0-8a3a-49f5-b49b-92e8e6d63507
function reaction(f, u, node, data)
    Tv=eltype(u)
    # buffer reactions
    rates = zeros(Tv, 5)
    ## in base
    ## CO2 + OH- <=> HCO3-
    rates[1] = kbf1 * u[ico2] * u[iohminus] - kbr1 * u[ihco3]
    ## HCO3- + OH- <=> CO3-- + H2O
    rates[2] = kbf2 * u[ihco3] * u[iohminus] - kbr2 * u[ico3]

    ## in acid
    ## CO2 + H20 <=> HCO3- + H+
    rates[3] = kaf1 * u[ico2] - kar1 * u[ihco3] * u[ihplus]
    ## HCO3- <=> CO3-- + H+
    rates[4] = kaf2 * u[ihco3] - kar2 * u[ico3] * u[ihplus]

    ## autoprotolyse
    rates[5] = kwf - kwr * u[ihplus] * u[iohminus]

    f[ihco3] -= rates[1] - rates[2] + rates[3] - rates[4]
    f[ico3] -= rates[2] + rates[4]
    f[ihplus] -= rates[3] + rates[4] + rates[5]
    f[iohminus] -= -rates[1] - rates[2] + rates[5]

    return nothing
end;

# ╔═╡ a7b7b5f9-148c-4361-95fb-00930f74c76b
md"""
### Boundary Conditions
"""

# ╔═╡ 7bddeb52-4a04-47a7-a848-3d5510bb074e
md"""
At the working electrode Robin boundary conditions are applied:

$-ε∇ϕ = C_{gap}~(ϕ_{we} - ϕ_{pzc} - ϕ^‡)$

where $ϕ^‡$ is the potential at the reaction plane.
"""

# ╔═╡ d45fa649-12b0-47cf-ae9c-a2baed705d8b
function halfcellbc(f, u, bnode, data)
    (; Γ_we, Γ_bulk, ϕ_we, iϕ) = data

    bulkbcondition(f, u, bnode, data; region = Γ_bulk)
    #boundary_dirichlet!(f, u, bnode; species = iϕ, region = Γ_we, value = ϕ_we)

    # Robin b.c. for the Poisson equation
    boundary_robin!(f, u, bnode, iϕ, Γ_we, C_gap, C_gap * (ϕ_we - ϕ_pzc))

    return nothing
end;

# ╔═╡ 73447419-de2a-47dd-ac38-e33823939f25
md"""
### Solver Control
"""

# ╔═╡ 704b905e-b456-4230-b9f7-573aad7660c8
solver_control = (
    max_round = 4,
    tol_round = 1.0e-8,
    reltol = 1.0e-8,
    abstol = 1.0e-9,
    verbose = "",
    maxiters = 20,
);

# ╔═╡ 2fbe5a81-063a-4e84-8799-c5acde5000cf
md"""
## Nernst-Planck Halfcell
"""

# ╔═╡ 115ead90-3694-4493-b1af-974965cba5db
grid = let
    hmin = 1.0e-6 * μm * 2.0^(-nref)
    hmax = 1.0 * μm * 2.0^(-nref)
    L = 80.0 * μm
    X = geomspace(0, L, hmin, hmax)
    simplexgrid(X)
end

# ╔═╡ e0f7b4c5-537e-474f-bb3b-3abd765dfdc8
celldata = ElectrolyteData(;
    nc = nc,
    na = 0,
    z = bulk.z,
    D = bulk.D,
    T = T,
    eneutral = false,
    κ = bulk.κ,
    c_bulk = bulk.c_bulk,
    Γ_we = 1,
    Γ_bulk = 2,
);

# ╔═╡ eff87b3e-0bc4-440f-93c9-9bd06e285a2e
@test isincompressible(celldata.c_bulk, celldata)

# ╔═╡ 4e31f737-1dc2-4505-8a1c-9ac523c220f3
@test iselectroneutral(celldata.c_bulk, celldata)

# ╔═╡ 5a7d3d0f-d337-44d2-b283-5b760bd3d14b
cell = PNPSystem(grid; bcondition = halfcellbc, reaction = reaction, celldata)

# ╔═╡ 40a4a273-e60c-40e8-8594-e041a95dcf5e
result = ivsweep(cell; voltages = vmin:vdelta:vmax, store_solutions = true, solver_control...)

# ╔═╡ e129c6b4-3e6a-4fbf-b6cc-5e87335142a4
md"""
### Results
"""

# ╔═╡ f73f98b6-2e0c-45a6-ab3c-611fa45c2d13
(~, default_index) = findmin(abs, result.voltages .+ 0.9 * ufac"V");

# ╔═╡ 9c609ace-f84d-45cb-8c65-6e29b548206d
md"""
$(@bind vindex PlutoUI.Slider(1:5:length(result.voltages), default=default_index))
"""

# ╔═╡ 287c0ac5-5f15-4f0d-9b9e-51ce873766e7
md"""
Potential at the working electrode 
$(vshow = result.voltages[vindex]; @sprintf("%+1.4f", vshow))
"""

# ╔═╡ 9636ac24-d78c-4f2b-b50b-eaa8313bca2a
md"""
Plot only the pH-value: $(@bind useonly_pH PlutoUI.CheckBox(default=false))
"""

# ╔═╡ 52051ed4-342e-48e2-a759-ee6c73378f3a
md"""
Compare with __Supplementary Figure 20__ in *Double layer charging driven carbon dioxide
adsorption limits the rate of electrochemical
carbon dioxide reduction on Gold* by __Ringe et al.__ in *Nature Communications*.
"""

# ╔═╡ 7b21cabc-651b-4bb0-be12-c0cccbc8b352
md"""
### Plotting Functions
"""

# ╔═╡ bf4d2747-0b36-41b5-a779-179f1897a4c6
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	function addplot(vis, sol, vshow)
		species = bulk.name
		colors = [:orange, :brown, :violet, :red, :blue, :green, :gray]
		
		scale = 1.0 / (mol / dm^3)
	    title = @sprintf("Φ_we=%+1.2f", vshow)
	
		if useonly_pH
			scalarplot!(vis, 
					    grid.components[XCoordinates] .+ 1.0e-14, 
					    log10.(sol[ihplus, :] * scale), 
					    color = colors[7],
					    label = species[7],
					    clear = true,
						title = title)
		else
			scalarplot!(vis, 
						grid.components[XCoordinates] .+ 1.0e-14, 
						log10.(sol[1, :] * scale), 
						color = colors[1],
						label = species[1],
						clear = true,
						title = title)
			for ia = 2:nc			
				scalarplot!(vis, 
						    grid.components[XCoordinates] .+ 1.0e-14, 
						    log10.(sol[ia, :] * scale), 
						    color = colors[ia],
						    label = species[ia],
						    clear = false,)
			end
		end
	end

	function plot1d(result, celldata, vshow)
		tsol 	= LiquidElectrolytes.voltages_solutions(result)
		vis 	= GridVisualizer(;
								 size 	= (600, 300),
								 clear 	= true,
								 legend 	= :rt,
								 limits 	= (-14, 2),
								 xlabel 	= "Distance from electrode [m]",
	 							 ylabel 	= "log c(aᵢ)", 
								 xscale 	= :log,)
	    addplot(vis, tsol(vshow), vshow)
		reveal(vis)
	end

	function plot1d(result, celldata)
    	tsol  	= LiquidElectrolytes.voltages_solutions(result)
		vis  	= GridVisualizer(; 
								 size 	= (600, 300),
								 clear 	= true,
							 	 legend = :rt,
								 limits = (-14, 2),
								 xlabel = "Distance from electrode [m]",
 								 ylabel = "log c(aᵢ)", 
								 xscale = :log,)
	
		vrange = result.voltages[1:5:end]
		movie(vis, file="bufferreactions.gif", framerate=6) do vis
		for vshow_it in vrange
			addplot(vis, tsol(vshow_it), vshow_it)
			reveal(vis)
		end
		end
		isdefined(Main, :PlutoRunner) && LocalResource("bufferreactions.gif")
	end
end
  ╠═╡ =#

# ╔═╡ 2856bb6b-7478-4696-98f4-ed5e233cb4ca
# ╠═╡ skip_as_script = true
#=╠═╡
plot1d(result, celldata, vshow)
  ╠═╡ =#

# ╔═╡ ceca65be-3614-4164-a71c-ca7d74e23bdc
# ╠═╡ skip_as_script = true
#=╠═╡
plot1d(result, celldata)
  ╠═╡ =#

# ╔═╡ 2a1ff350-0da1-4fa0-9e8d-2975c6d2d81e
TableOfContents(title = "📚 Table of Contents", indent = true, depth = 4, aside = true)

# ╔═╡ 4f56ddcf-5001-4f0f-b7de-b9c311f749a4
begin
    hrule() = html"""<hr>"""
    function highlight(mdstring, color)
        return htl"""<blockquote style="padding: 10px; background-color: $(color);">$(mdstring)</blockquote>"""
    end

    macro important_str(s)
        return :(highlight(Markdown.parse($s), "#ffcccc"))
    end
    macro definition_str(s)
        return :(highlight(Markdown.parse($s), "#ccccff"))
    end
    macro statement_str(s)
        return :(highlight(Markdown.parse($s), "#ccffcc"))
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

# ╔═╡ c875fc12-4c92-42f9-82ef-34ba242899c9
hrule()


# ╔═╡ Cell order:
# ╠═dbaff000-beeb-46c3-bbab-7c7e3a3e9791
# ╠═3b5a64e6-24cd-423e-aad8-8f400b338867
# ╟─eea5aa29-85d4-4644-bfe4-4b5244eeb045
# ╟─8d717a52-d9d5-491b-8ce8-98da978bbca5
# ╠═b263c6bf-3606-4c3d-9fe1-5d074b294519
# ╟─4c355a5e-ccb8-4d23-be5c-adc1aa5f8c96
# ╟─327ed52b-6c81-4c2f-8b79-824017351020
# ╠═e8de40fe-472e-437f-92c8-5f75f5d58601
# ╟─28309673-e34e-4bf3-9c53-d04a283364ef
# ╠═2a787207-55ab-4e5b-a9f2-22261ee28a42
# ╠═960fe09f-8bfe-4773-a57e-cbd7774050ce
# ╟─e969bc7f-f0f3-4f70-a522-133ee56cce5e
# ╟─3161211b-3188-4e4f-a12e-5a16cfd5e932
# ╟─18e13140-2693-4a2b-b8e2-24ae437414cd
# ╟─c1892a61-dd26-4173-8675-8136c8360469
# ╠═12c99cb0-8a3a-49f5-b49b-92e8e6d63507
# ╟─a7b7b5f9-148c-4361-95fb-00930f74c76b
# ╟─7bddeb52-4a04-47a7-a848-3d5510bb074e
# ╠═d45fa649-12b0-47cf-ae9c-a2baed705d8b
# ╟─73447419-de2a-47dd-ac38-e33823939f25
# ╠═704b905e-b456-4230-b9f7-573aad7660c8
# ╟─2fbe5a81-063a-4e84-8799-c5acde5000cf
# ╠═115ead90-3694-4493-b1af-974965cba5db
# ╠═e0f7b4c5-537e-474f-bb3b-3abd765dfdc8
# ╠═eff87b3e-0bc4-440f-93c9-9bd06e285a2e
# ╠═4e31f737-1dc2-4505-8a1c-9ac523c220f3
# ╠═5a7d3d0f-d337-44d2-b283-5b760bd3d14b
# ╠═40a4a273-e60c-40e8-8594-e041a95dcf5e
# ╟─e129c6b4-3e6a-4fbf-b6cc-5e87335142a4
# ╟─f73f98b6-2e0c-45a6-ab3c-611fa45c2d13
# ╟─287c0ac5-5f15-4f0d-9b9e-51ce873766e7
# ╟─9c609ace-f84d-45cb-8c65-6e29b548206d
# ╟─9636ac24-d78c-4f2b-b50b-eaa8313bca2a
# ╟─52051ed4-342e-48e2-a759-ee6c73378f3a
# ╠═2856bb6b-7478-4696-98f4-ed5e233cb4ca
# ╠═ceca65be-3614-4164-a71c-ca7d74e23bdc
# ╟─7b21cabc-651b-4bb0-be12-c0cccbc8b352
# ╠═bf4d2747-0b36-41b5-a779-179f1897a4c6
# ╟─c875fc12-4c92-42f9-82ef-34ba242899c9
# ╠═2a1ff350-0da1-4fa0-9e8d-2975c6d2d81e
# ╟─4f56ddcf-5001-4f0f-b7de-b9c311f749a4
