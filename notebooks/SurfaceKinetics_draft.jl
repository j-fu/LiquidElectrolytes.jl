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

# ╔═╡ 5874f7fa-1bf8-11ee-0ace-859064fd872c
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "docs")) # hide
    using Revise
    using LessUnitful
    using ExtendableGrids,GridVisualize
    using LiquidElectrolytes
    using VoronoiFVM
    using InteractiveUtils
    using ForwardDiff
    using PlutoUI, HypertextLiteral
	using DataFrames
	using Test
	using Printf
	using Interpolations
    if isdefined(Main,:PlutoRunner)
		using CairoMakie,Colors 
   		default_plotter!(CairoMakie)
		CairoMakie.activate!(type="png")
    end
end

# ╔═╡ 39e030bf-2280-4c24-9f7f-3f4c0b1ca2b0
md"""
# Surface kinetics
"""

# ╔═╡ 270c17fd-a168-4b74-b9ca-362a647e6e89
md"""
A toy model for a electrochemical half cell including a simplified microkinetic model for the electrochemical surface reactions is simulated. 

In an electronetural solvent (e.g. H₂O) a salt A⁺B⁻ fully dissociates into solvated A⁺-ions and B⁻-ions. At the working electrode A⁺-ions are reduced:

$A^+_{aq} + e^-\to A_{aq}$

and A-atoms are oxidized

$A_{aq} \to e^- + A^+_{aq}$.

A simplified microkinetic model for the surface reactions is applied that contains the following steps:

1. Adsorption: $A^+ \rightleftharpoons A^+_{ad}$
2. Electron transfer: $A^+_{ads} + e^- \rightleftharpoons A_{ads}$
3. Desorption: $A_{ads} \rightleftharpoons A_{aq}$.

A mean-field approach using the law of mass action is applied. No transition states between the reactants are assumed. The reaction rate constants are model by the Arrhenius relation with pre-exponential factors of size 10¹³. A linear scaling of the adsorption energies with the (excess) surface charge at the electrode is assumed:

$ΔG_{ads}(σ) = ΔG_{ads}(σ = 0) + b * σ$

where the model $\sigma = -\varepsilon \nabla \phi|_{\text{we}}\cdot \vec{n}|_{\text{we}} = C_{\text{gap}} (ϕ_{\text{we}} - \phi_{\text{pzc}} - \phi^\ddagger)$ (with $\phi^\ddagger$ being the potential at the reaction plane) is used.

The potential dependence of the electron transfer reaction is given by:

$ΔG_{\text{rxn}}(U) = ΔG_{\text{rxn}}(U = 0) + ϕ_{we}$
"""

# ╔═╡ f73b2304-0c5a-45b4-9807-1fe068ef4dc9
pkgdir(LiquidElectrolytes)

# ╔═╡ b2d028c5-3d44-4ac5-ae11-cd0416d89245
md"""
## Setup
"""

# ╔═╡ de330a48-dc23-480a-ad17-07b1e7dae042
md"""
### Units
"""

# ╔═╡ c0f31ee2-fb00-4bd0-9762-81d3791695de
begin
	@unitfactors eV V μF cm K mol m A s dm μm nm;
	@phconstants N_A e R k_B
	F = N_A * e
end;

# ╔═╡ a876c117-2d55-4021-83d6-8c05aaafd579
md"""
### Data
"""

# ╔═╡ 2e585c01-717b-43c0-a79b-2c574c6fdd11
begin
	const T = (273.15 + 25) * K
	const C_gap = 20 * μF / cm^2
	const ϕ_pzc = 0.1 * V

	const S = 1.0e-5 / N_A * (1.0e10)^2 * mol / m^2

	const ΔG_ads_σ0_aplus = -0.01 * eV
	const ΔG_ads_σ0_a = -0.005 * eV

	const b_ads_aplus = 0.5 * eV * m^2 / (A * s)
	const b_ads_a = 0.25 * eV * m^2 / (A * s)

	const ΔG_rxn_U0 = -0.1 * eV
	
	const vmin = -0.2 * V
	const vmax =  0.3 * V
	const vdelta = 0.1 * V

	const nref = 0

	const scheme = :μex

	const iaplus      = 1
    const ibminus     = 2
    const ia          = 3
    const iaplus_ads  = 4
    const ia_ads      = 5
end;

# ╔═╡ 7293cb64-b3a4-4ec3-aff5-b62104716e1b
const bulk = DataFrame(
  :name 	=> ["A⁺", "B⁻", "A"],
  :z 		=> [1, -1, 0],
  :D 		=> [2.0e-9, 2.0e-9, 2.0e-9] * m^2/s,
  :κ    	=> [8.0, 4.0, 2.0],
  :c_bulk 	=> [0.1, 0.1, 0.1] * mol/dm^3
)

# ╔═╡ 32704bd9-840c-48df-9174-ccec27dea865
md"""
### Surface Reactions
"""

# ╔═╡ 54e39b03-e646-474a-8ba2-99979a99b201
function rateconstants!(k, σ, ϕ_we, ϕ)
	(kf, kr) = k

	#println("σ 	   (at ϕ_we = $(@sprintf("%+1.2f", ϕ_we))) = $(ForwardDiff.value(σ))")
	# A+_aq <-> A+_ads,	                    #1
	ΔG_ads_aplus = ΔG_ads_σ0_aplus + b_ads_aplus * σ
	#println("ΔG_ads (at ϕ_we = $(@sprintf("%+1.2f", ϕ_we))) = $(ForwardDiff.value(ΔG_ads_aplus))")

	kf[1] = 1.0e13 * exp(-max(ΔG_ads_aplus, 0.0) / (k_B * T))
	kr[1] = 1.0e13 * exp(-max(-ΔG_ads_aplus, 0.0) / (k_B * T))

	# A+_ads + e- <-> A_ads,                #2          
	ΔG_rxn = ΔG_rxn_U0 + (ϕ_we - ϕ) * eV
	#println("ΔG_rxn (at ϕ_we = $(@sprintf("%+1.2f", ϕ_we))) = $(ForwardDiff.value(ΔG_rxn))")

	kf[2] = 1.0e13 * exp(-max(ΔG_rxn, 0.0) / (k_B * T))
	kr[2] = 1.0e13 * exp(-max(-ΔG_rxn, 0.0) / (k_B * T))

	# A_ads <-> A_aq                        #3
	ΔG_ads_a = ΔG_ads_σ0_a + b_ads_a * σ

	kf[3] = 1.0e13 * exp(-max(ΔG_ads_a, 0.0) / (k_B * T))
	kr[3] = 1.0e13 * exp(-max(-ΔG_ads_a, 0.0) / (k_B * T))
end;

# ╔═╡ bf7e4757-0ac3-47db-a518-68db32217800
function breactions(
	f,
	u::VoronoiFVM.BNodeUnknowns{Tval, Tv, Tc, Tp, Ti}, 
	bnode,
	data
) where {Tval, Tv, Tc, Tp, Ti}
	(; Γ_we, ϕ_we, iϕ) = data
	
	if bnode.region == Γ_we
		
		σ = C_gap * (ϕ_we - u[iϕ] - ϕ_pzc)         

		kf = zeros(Tval, 3)
		kr = zeros(Tval, 3)

		rateconstants!((kf, kr), σ, ϕ_we, u[iϕ])

		rates = zeros(Tval, 3)

		θ_free = 1 - u[iaplus_ads]  - u[ia_ads]
		rates[1] = kf[1] * u[iaplus] * θ_free - kr[1] * u[iaplus_ads]
		rates[2] = kf[2] * u[iaplus_ads] - kr[2] * u[ia_ads]
		rates[3] = kf[3] * u[ia_ads] - kr[3] * u[ia] * θ_free

		# bulk species
		f[iaplus] += -rates[1] * S
		f[ia] += rates[3] * S
		
		# surface species
		f[iaplus_ads] += rates[1] - rates[2]
		f[ia_ads] += rates[2] - rates[3]
	end
	nothing
end;

# ╔═╡ c32bf149-320b-47a3-9e1b-fd6a29d39834
function halfcellbc(f,u, bnode,data)
	(; Γ_we, Γ_bulk, ϕ_we, iϕ) = data

	bulkbcondition(f, u, bnode, data; region = Γ_bulk)
	boundary_robin!(f, u, bnode, iϕ, Γ_we , C_gap, C_gap * (ϕ_we - ϕ_pzc))
	breactions(f, u, bnode, data)
	
	nothing
end;

# ╔═╡ abc3b243-c6af-4bd8-a2fe-5ac488a5cd3e
md"""
### Solver Control
"""

# ╔═╡ ba2cf338-ec72-4129-829d-698a8473787a
solver_control = (max_round = 4,
                  tol_round = 1.0e-8,
                  reltol 	= 1.0e-8,
                  abstol 	= 1.0e-9,
                  verbose 	= "",
                  maxiters 	= 20);

# ╔═╡ 40db90b2-5ee3-4773-9d0b-bdaa1c095cf2
md"""
## Nernst-Planck Half-Cell
"""

# ╔═╡ a12023d9-44b6-4db9-9b37-05b3c371327b
grid = let
	hmin    = 1.0e-2 * nm * 2.0^(-nref)
    hmax    = 1.0 * nm * 2.0^(-nref)
    L       = 20.0 * nm
    X       = geomspace(0, L, hmin, hmax)
    simplexgrid(X)
end

# ╔═╡ 50e5aebb-e42e-4fe3-8f39-601b9532c63f
celldata = ElectrolyteData(; nc 	= 3,
							 na 	= 2,
							 z 		= bulk.z,
							 D 		= bulk.D,
							 T 		= T,
							 eneutral = false,
							 κ 		= bulk.κ,
							 c_bulk = bulk.c_bulk,
							 Γ_we 	= 1,
							 Γ_bulk = 2,
							 scheme);

# ╔═╡ 0824a2e2-0753-4891-b265-d919b694df31
@test isincompressible(celldata.c_bulk, celldata)

# ╔═╡ 9f2397ad-2842-459b-b584-5c0595ff9f02
@test iselectroneutral(celldata.c_bulk, celldata)

# ╔═╡ 5b5be22a-f713-4cf5-bc93-77864d198149
cell = PNPSystem(grid; bcondition=halfcellbc, celldata);

# ╔═╡ e3248c10-e721-4c76-997f-34bd4370591c
md"""
## Results
"""

# ╔═╡ 6fe6bb3a-4182-4a70-9dcc-947584e0a83e
result = ivsweep(cell; voltages=vmin:vdelta:vmax, store_solutions=true, solver_control...);

# ╔═╡ 5a59343c-6f2a-4ffd-935a-f265ce10c70c
@bind vshow PlutoUI.Slider(range(extrema(result.voltages)...; length = 101),
                           show_value = true) # hide

# ╔═╡ d1fec937-c3fa-41a4-bfbd-f73bcc2e6015
md"""
### Plotting functions
"""

# ╔═╡ f53c55dd-cb8f-4bb4-91f2-3311abe1a9d3
curr(J, ix) = [F * j[ix] for j in J]

# ╔═╡ c3385ba7-357a-41bb-87d3-dc34739669c5
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
                curr(result.j_bulk, iaplus);
                linestyle = :dash,
                label = "A⁺, bulk",
                color = :red)
    scalarplot!(vis,
                volts,
                curr(result.j_we, iaplus);
                color = :red,
                clear = false,
                linestyle = :solid,
                label = "A⁺, we")

    scalarplot!(vis,
                volts,
                curr(result.j_bulk, ia);
                linestyle = :dash,
                label = "A, bulk",
                color = :green,
                clear = false)
    scalarplot!(vis,
                volts,
                curr(result.j_we, ia);
                color = :green,
                clear = false,
                linestyle = :solid,
                label = "A, we")

    scalarplot!(vis,
                volts,
                curr(result.j_bulk, ibminus);
                linestyle = :dash,
                label = "B⁻, bulk",
                color = :green,
                clear = false)
    scalarplot!(vis,
                volts,
                curr(result.j_we, ibminus);
                color = :green,
                clear = false,
                linestyle = :solid,
                label = "B⁻, we")

    reveal(vis)
end

# ╔═╡ a83b42f5-c7a8-48c1-af96-c928c80b8f8c
plotcurr(result)

# ╔═╡ 0b31ab30-16f4-4c45-abe9-792fa12400b4
begin
	function addplot(vis, sol, vshow, curr)
		scale = 1.0 / (mol / dm^3)
	    title = @sprintf("ϕ_we=%+1.2f, I=%+1.2f", vshow, curr)
	    c0 = solventconcentration(sol, celldata)
	    
		ϵ = 1.0e-14
		scalarplot!(vis, 
					grid.components[XCoordinates] .+ ϵ, 
					sol[ia, :] * scale; 
					color = :green, 
					label = "A",
					clear = true,
					title)
	    scalarplot!(vis,
	                grid.components[XCoordinates] .+ ϵ,
	                sol[ibminus, :] * scale;
	                color = :gray,
	                clear = false,
	                label = "B⁻")
	    scalarplot!(vis,
	                grid.components[XCoordinates] .+ ϵ,
	                sol[iaplus, :] * scale;
	                color = :red,
	                clear = false,
	                label = "A⁺")
	    scalarplot!(vis, 
					grid.components[XCoordinates] .+ ϵ, 
					c0 * scale; 
					color = :blue, 
					clear = false,
	                label = "H2O")
	end

	function plot1d(result, celldata, vshow)
		tsol   = LiquidElectrolytes.voltages_solutions(result)
		vinter = linear_interpolation(result.voltages, [j[ia] for j in result.j_we])
		
		vis = GridVisualizer(;
					 size = (600, 250),
					 yscale = :log,
					 limits = (1.0e-5, 1.0e2),
					 xlimits = (1.0e-11, 20 * nm),
					 legend = :rt)
		addplot(vis, tsol(vshow), vshow, vinter(vshow))
		reveal(vis)
	end
	
	function plot1d(result, celldata)
		tsol   = LiquidElectrolytes.voltages_solutions(result)
		vinter = linear_interpolation(result.voltages, [j[ia] for j in result.j_we])
		vrange = result.voltages[1:1:end]
		
		vis = GridVisualizer(;
					 size = (600, 250),
					 yscale = :log,
					 limits = (1.0e-5, 1.0e2),
					 xlimits = (1.0e-11, 20 * nm),
					 legend = :rt)
		movie(vis, file="surfacereactions.gif", framerate=6) do vis
		for vshow_it in vrange
			addplot(vis, tsol(vshow_it), vshow_it, vinter(vshow_it))
			reveal(vis)
		end
		end
		isdefined(Main, :PlutoRunner) && LocalResource("surfacereactions.gif")
	end
end

# ╔═╡ 7002cb8e-8010-4b60-a707-7537c494636f
plot1d(result, celldata, vshow)

# ╔═╡ dbc76012-4458-4ef2-897a-dd8e4ced4931
plot1d(result, celldata)

# ╔═╡ 8812965c-9ab9-47c7-9ca2-7d08c8157a69
function cplot(cell, result)
    scale = 1.0 / (mol / dm^3)
    tsol=LiquidElectrolytes.voltages_solutions(result)
    j_we=result.j_we
    currs = curr(j_we, ia)
    vis = GridVisualizer(; resolution = (1200, 400), layout = (1, 3),
                         gridscale = 1.0e9)
    xmax = 10 * nm
    xlimits = [0, xmax]
    aspect = 2 * xmax / (tsol.t[end] - tsol.t[begin])
    scalarplot!(vis[1, 1],
                cell,
                tsol;
                scale,
                species = ia,
                aspect,
                xlimits,
                title = "A",
                colormap = :summer)
    scalarplot!(vis[1, 2],
                cell,
                tsol;
                species = iaplus,
                aspect,
                scale,
                xlimits,
                title = "A⁺",
                colormap = :summer)
    scalarplot!(vis[1, 3],
                tsol[ia, 1, :] * scale,
                tsol.t;
                label = "A",
                xlabel = "c",
                color = :green,
                clear = false)
    scalarplot!(vis[1, 3],
                tsol[iaplus, 1, :] * scale,
                tsol.t;
                title = "c(0)",
                xlabel = "c",
                ylabel = "V",
                label = "A⁺",
                color = :red,
                clear = false)
    scalarplot!(vis[1, 3],
                tsol[ibminus, 1, :] * scale,
                tsol.t;
                label = "B⁻",
                color = :blue,
                clear = false,
                legend = :ct)
    reveal(vis)
end

# ╔═╡ fb0530a9-9ecf-4eda-83c1-89de6ee93d0d
cplot(cell, result)

# ╔═╡ 855fbd2c-5dcd-47a7-84b8-c859b3f80ba8
TableOfContents()

# ╔═╡ 2ca1204a-56a6-4d9a-bc8a-46145382430c
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

# ╔═╡ f1cfba83-bea2-43ae-9dc8-fd2788a4b937
hrule()


# ╔═╡ Cell order:
# ╟─39e030bf-2280-4c24-9f7f-3f4c0b1ca2b0
# ╟─270c17fd-a168-4b74-b9ca-362a647e6e89
# ╠═5874f7fa-1bf8-11ee-0ace-859064fd872c
# ╠═f73b2304-0c5a-45b4-9807-1fe068ef4dc9
# ╟─b2d028c5-3d44-4ac5-ae11-cd0416d89245
# ╟─de330a48-dc23-480a-ad17-07b1e7dae042
# ╠═c0f31ee2-fb00-4bd0-9762-81d3791695de
# ╟─a876c117-2d55-4021-83d6-8c05aaafd579
# ╠═2e585c01-717b-43c0-a79b-2c574c6fdd11
# ╠═7293cb64-b3a4-4ec3-aff5-b62104716e1b
# ╟─32704bd9-840c-48df-9174-ccec27dea865
# ╠═c32bf149-320b-47a3-9e1b-fd6a29d39834
# ╠═bf7e4757-0ac3-47db-a518-68db32217800
# ╠═54e39b03-e646-474a-8ba2-99979a99b201
# ╟─abc3b243-c6af-4bd8-a2fe-5ac488a5cd3e
# ╠═ba2cf338-ec72-4129-829d-698a8473787a
# ╟─40db90b2-5ee3-4773-9d0b-bdaa1c095cf2
# ╠═a12023d9-44b6-4db9-9b37-05b3c371327b
# ╠═50e5aebb-e42e-4fe3-8f39-601b9532c63f
# ╠═0824a2e2-0753-4891-b265-d919b694df31
# ╠═9f2397ad-2842-459b-b584-5c0595ff9f02
# ╠═5b5be22a-f713-4cf5-bc93-77864d198149
# ╟─e3248c10-e721-4c76-997f-34bd4370591c
# ╠═6fe6bb3a-4182-4a70-9dcc-947584e0a83e
# ╠═a83b42f5-c7a8-48c1-af96-c928c80b8f8c
# ╟─5a59343c-6f2a-4ffd-935a-f265ce10c70c
# ╠═7002cb8e-8010-4b60-a707-7537c494636f
# ╠═dbc76012-4458-4ef2-897a-dd8e4ced4931
# ╠═fb0530a9-9ecf-4eda-83c1-89de6ee93d0d
# ╟─d1fec937-c3fa-41a4-bfbd-f73bcc2e6015
# ╟─f53c55dd-cb8f-4bb4-91f2-3311abe1a9d3
# ╟─c3385ba7-357a-41bb-87d3-dc34739669c5
# ╟─0b31ab30-16f4-4c45-abe9-792fa12400b4
# ╟─8812965c-9ab9-47c7-9ca2-7d08c8157a69
# ╟─f1cfba83-bea2-43ae-9dc8-fd2788a4b937
# ╟─855fbd2c-5dcd-47a7-84b8-c859b3f80ba8
# ╟─2ca1204a-56a6-4d9a-bc8a-46145382430c
