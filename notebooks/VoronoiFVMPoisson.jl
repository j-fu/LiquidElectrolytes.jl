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

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
begin

    import Pkg # hide
    Pkg.activate(joinpath(@__DIR__, "..", "docs")) # hide
    using Revise # hide

    using PlutoUI, HypertextLiteral, UUIDs
    using LinearAlgebra
    using Interpolations
    using VoronoiFVM, GridVisualize, ExtendableGrids
    using LaTeXStrings
    using LessUnitful, Unitful

    if isdefined(Main, :PlutoRunner)
        using CairoMakie
        default_plotter!(CairoMakie)
        CairoMakie.activate!(type = "svg")
    end

end

# ╔═╡ 66dc2136-a601-4291-9aa5-81616a4efdd7
md"""
# VoronoiFVMPoisson.jl
"""

# ╔═╡ c8ac41dd-df4d-47f6-909e-4e18e65ca994
md"""
Realization of Poisson-Boltzmann and Poisson-Bikerman equations based on plain VoronoiFVM.jl (without LiquidElectrolytes). This notebook is meant for demonstation purposes and may be moved to VoronoiFVM.jl.
"""

# ╔═╡ b181a82b-46e7-4dae-b0b1-6886dd40886a
md"""
## Constants, Units & Parameters
"""

# ╔═╡ 5a5e3586-2009-4381-ae9d-a0bd34539b65
md"""
See LessUnitful.jl  `ph` and `ufac` and Unitful.jl for `u`
"""

# ╔═╡ b6338637-1f3c-4410-8e82-3afdaf603656
begin
    const ε = 78.49
    const ε_0 = ph"ε_0"
    const F = ph"N_A" * ph"e"
    const K = ufac"K"
    const nm = ufac"nm"
    const dm = ufac"dm"
    const V = ufac"V"
    const mol = ufac"mol"
end;

# ╔═╡ 371df98d-9f55-4244-8243-e6b597cc8d84
md"""
Parameters:
"""

# ╔═╡ 10be79ba-9ab9-4e6f-8568-2de93ea77964
begin
    const molarity = 0.01
    const c_bulk = [molarity, molarity] * mol / dm^3 # bulk ion molar concentrations
    const z = [-1, 1] # species charge numbers
    const c̄ = 55.508mol / dm^3 # summary molar concentration
    const nref = 0 # grid refinement level
    const T = 298.15K # temperature
    const L = 20.0nm # computational domain size
    const vmax = 1.0 * V # maximum voltage
end;

# ╔═╡ 834c81c8-9035-449f-b331-73bbee87756c
md"""
Check for bulk electroneutrality:
"""

# ╔═╡ 7f209979-0776-40ab-9b2b-3b2d0145fa2c
@assert dot(c_bulk, z) == 0

# ╔═╡ 80a42b07-d6b1-41d6-bf1d-66a01eeb5019
md"""
Derived data:
"""

# ╔═╡ fb653542-7084-4e6f-87b9-9ab29e89ba04
begin
    const RT = ph"R" * T
    const c0_bulk = c̄ - sum(c_bulk) # solvent bulk molar concentration
    const hmin = 1.0e-1 * nm * 2.0^(-nref) # grid size at working electrode
    const hmax = 1.0 * nm * 2.0^(-nref) # grid size at bulk
    const l_debye = sqrt(ε * ε_0 * RT / (F^2 * c_bulk[1])) # Debye length
    const dlcap0 = sqrt(2 * ε * ε_0 * F^2 * c_bulk[1] / RT) # Double layer capacitance at point of zero charge (0V)
end;

# ╔═╡ c57e82bb-79d8-4553-b31d-d1448c38c649
md"""
Debye length= $(l_debye |> x->round(x,sigdigits=5) |>u"nm")
"""

# ╔═╡ e2179a6c-cc1c-4850-89d4-d3f87fd5e6ee
md"""
Double layer capacitance at zero voltage for symmetric binary electrolyte = 
$(dlcap0 |> x->round(x,sigdigits=5) |>u"μF/cm^2")
"""

# ╔═╡ 7da889ce-9c6b-4abc-b19d-6311aacc32b1
X = geomspace(0, L, hmin, hmax)

# ╔═╡ a2a8132a-d54d-48a7-aa02-ff83593f0e16
grid = simplexgrid(X)

# ╔═╡ db023dac-5ce0-4126-bfd0-9036ccec11f4
gridplot(grid, size = (600, 200))

# ╔═╡ 7ee7be9f-038d-4080-af7c-01ae8f7d2682
md"""
## Nonlinear Poisson equation
"""

# ╔═╡ 40282ce8-e032-4524-9a20-971826a6bb97
md"""
Let ``\bar c`` be the summary concentration, and ``y_i`` be the solute molar fractions, such that for the species molar concentrations one has ``c_i=\bar c y_i``. Then the Poisson equation states:
```math
\begin{aligned}
  -\nabla \cdot εε_0 \nabla \phi &= F\sum_{i=1}^n z_i \bar c y_i\\
\end{aligned}
```
"""

# ╔═╡ 1b13b664-9d5e-43a1-a7b2-2a1c6da4c947
md"""
Un VoronoiFVM, the left hand side is expressed via:
"""

# ╔═╡ 5befbcff-2e7c-48db-a327-27dbb6779d63
function pb_flux(y, u, edge, data)
    return y[1] = ε * ε_0 * (u[1, 1] - u[1, 2])
end

# ╔═╡ 91ecd9f8-5823-4ba0-9e39-a644f24424ef
md"""
It remains to state the dependecy of the mole fractions on the potential.
"""

# ╔═╡ f8914d81-117d-4d9e-8f89-27a3045c8e7f
md"""
We choose Dirichlet boundary conditions at the electrode (index 1) and the bulk (index 2):
"""

# ╔═╡ 2a7c5a4f-48b0-4995-a37b-70fd3d5bf853
const bcvals = [0.0, 0.0]

# ╔═╡ ddbb2a10-f22e-4474-b02d-e7b514852e58
function pb_bc(y, u, bnode, data)
    boundary_dirichlet!(y, u, bnode, species = 1, region = 2, value = bcvals[2])
    return boundary_dirichlet!(y, u, bnode, species = 1, region = 1, value = bcvals[1])
end

# ╔═╡ a4a4c6a4-ea90-4d23-b13a-2020790b2889
md"""
## Classical Poisson-Boltzmann

```math
\begin{aligned}
   y_i&=  \frac{c_{i,bulk}}{\bar c}\exp\left(z_i\phi \frac{F}{RT}\right)
\end{aligned}
```
"""

# ╔═╡ 286174fc-9a64-426a-9d3c-a96b39f7ea0c
pbo_molfrac(ϕ, i) = (c_bulk[i] / c̄) * exp(z[i] * ϕ * F / RT)

# ╔═╡ 39996246-cb69-4c59-802a-d00804d84b0b
function pbo_spacecharge(ϕ)
    sumyz = zero(ϕ)
    for i in 1:length(z)
        sumyz += z[i] * pbo_molfrac(ϕ, i)
    end
    return F * c̄ * sumyz
end

# ╔═╡ 95e49176-6c34-4c89-b156-5972080f82bb
function pbo_reaction(y, u, node, data)
    return y[1] = pbo_spacecharge(u[1])
end

# ╔═╡ 91a661ee-9c0c-496e-b6fb-4bb692cb7255
pbo_system = VoronoiFVM.System(
    grid;
    reaction = pbo_reaction,
    flux = pb_flux,
    bcondition = pb_bc,
    species = [1]
)

# ╔═╡ 0c77e9c2-1242-4694-bacf-3b6c9856997a
function sweep(system; nsweep = 100, δ = 1.0e-4)
    voltages = range(0, vmax, length = nsweep + 1)
    solutions = []
    dlcaps = zeros(0)
    for v in voltages
        bcvals[1] = v
        if length(solutions) == 0
            inival = 0.0
        else
            inival = solutions[end]
        end
        push!(solutions, solve(system; inival))
        q = integrate(system, system.physics.reaction, solutions[end])
        bcvals[1] = v + δ
        solδ = solve(system; inival = solutions[end])
        qδ = integrate(system, system.physics.reaction, solδ)
        push!(dlcaps, (qδ[1] - q[1]) / δ)
    end
    return voltages, dlcaps, TransientSolution(solutions, voltages)
end

# ╔═╡ 163c5395-726f-4874-9304-ab9a46cc8bf9
pbo_volts, pbo_caps, pbo_sols = sweep(pbo_system)

# ╔═╡ db5b6820-5a53-465c-b380-66756fd722a6
md"""
## "Poisson-Bikerman"
"""

# ╔═╡ 11bbb453-0a80-4ec4-b5a9-c044ebce8146
md"""
Here we use the  expression
```math
\begin{aligned}
   y_i&= \frac{\frac{c_{i,bulk}}{c_{0,bulk}}\exp\left(z_i\phi \frac{F}{RT}\right)}{1+\sum_{j=1}^n \frac{c_{j,bulk}}{c_{0,bulk}}\exp\left(z_i\phi \frac{F}{RT}\right)}
\end{aligned}
```
"""

# ╔═╡ c2500d41-53e8-4ac7-8e52-5107cf0335b2
function pbi_molfrac(ϕ, j)
    N = length(z)
    denom = one(ϕ)
    for i in 1:length(z)
        denom += exp(z[i] * ϕ * F / RT) * c_bulk[i] / c0_bulk
    end
    return exp(z[j] * ϕ * F / RT) * c_bulk[j] / (c0_bulk * denom)
end

# ╔═╡ cb8b5a32-c414-4768-b250-2398cc432881
function pbi_spacecharge(ϕ)
    sumyz = zero(ϕ)
    for i in 1:length(z)
        sumyz += z[i] * pbi_molfrac(ϕ, i)
    end
    return F * c̄ * sumyz
end

# ╔═╡ 77e92a09-66e8-4502-8115-c892eef1f597
function pbi_reaction(y, u, node, data)
    return y[1] = pbi_spacecharge(u[1])
end

# ╔═╡ 0b628db9-1719-4e1a-b87b-92d3e2808e5c
pbi_system = VoronoiFVM.System(grid; reaction = pbi_reaction, flux = pb_flux, bcondition = pb_bc, species = [1])

# ╔═╡ db911ce1-c0f1-4261-aee7-54156ed03248
pbi_volts, pbi_caps, pbi_sols = sweep(pbi_system)

# ╔═╡ e373824d-6e0d-44c8-8c8a-192cccf10298
function plotcdl(pbvolts, pbcaps, pbv)
    vc = linear_interpolation(pbvolts, pbcaps)
    cdl = round(vc(pbv) / ufac"μF/cm^2", sigdigits = 4)
    vis = GridVisualizer(size = (500, 200), xlabel = "Δϕ", ylabel = L"C_{dl}/(μF/cm^2)", legend = :lt)
    scalarplot!(vis, pbvolts, pbcaps / ufac"μF/cm^2")
    scalarplot!(
        vis, [pbv], [cdl], markershape = :circle, clear = false, markersize = 20,
        label = "$(cdl)"
    )
    return reveal(vis)
end

# ╔═╡ a290a533-9639-465f-b6ef-66c28d6136e8
function plotsols(pbsols, pbv, pb_molfrac)
    vis = GridVisualizer(size = (700, 200), layout = (1, 2), xlabel = "x/nm")
    ϕ = pbsols(pbv)[1, :]
    cp = c̄ * pb_molfrac.(ϕ, 1) / ufac"mol/dm^3"
    cm = c̄ * pb_molfrac.(ϕ, 2) / ufac"mol/dm^3"
    scalarplot!(vis[1, 1], X / nm, ϕ, color = :green, ylabel = "ϕ/V")
    scalarplot!(vis[1, 2], X / nm, cp, color = :red, ylabel = "c/(mol/L)")
    scalarplot!(vis[1, 2], X / nm, cm, color = :blue, clear = false)
    return reveal(vis)
end

# ╔═╡ f9b4d4dc-7def-409f-b40a-f4eba1163741
TableOfContents()

# ╔═╡ 7a93e9a8-8a2d-4b11-84ef-691706c0eb0f
begin
    hrule() = html"""<hr>"""
    highlight(mdstring, color) =
    htl"""<blockquote style="padding: 10px; background-color: $(color);">$(mdstring)</blockquote>"""

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
        	/* Headers */
             h1{background-color:#dddddd;  padding: 10px;}
             h2{background-color:#e7e7e7;  padding: 10px;}
             h3{background-color:#eeeeee;  padding: 10px;}
             h4{background-color:#f7f7f7;  padding: 10px;}

    		/* "Terminal"  */
    	     pluto-log-dot-sizer  { max-width: 655px;}
             pluto-log-dot.Stdout { background: #002000;
    	                            color: #10f080;
                                    border: 6px solid #b7b7b7;
                                    min-width: 18em;
                                    max-height: 300px;
                                    width: 675px;
                                    overflow: auto;
     	                           }
            /* Standard cell width etc*/
    		main {
       			flex: 1;
    		    max-width: calc(700px + 25px + 6px); /* 700px + both paddings */
        		padding-top: 0px;
        		padding-bottom: 4rem;
        		padding-left: 25px;
        		padding-right: 6px;
        		align-content: center;
        		width: 100%;
    			}

           /* Cell width for slides*/
    		xmain {
    			margin: 0 auto;
    			max-width: 750px;
        		padding-left: max(20px, 3%);
        		padding-right: max(20px, 3%);
    	        }

    	
    	
        </style>
    """
end

# ╔═╡ 9b5f389f-b105-4610-bab7-f79305fedc31
hrule()

# ╔═╡ 5beb3a0d-e57a-4aea-b7a0-59b8ce9ff5ce
hrule()

# ╔═╡ 86e1fffe-cf70-4d93-8131-153f0e98c94b
begin
    function floataside(text::Markdown.MD; top = 1)
        uuid = uuid1()
        return @htl(
            """
            		<style>


            		@media (min-width: calc(700px + 30px + 300px)) {
            			aside.plutoui-aside-wrapper-$(uuid) {

            	color: var(--pluto-output-color);
            	position:fixed;
            	right: 1rem;
            	top: $(top)px;
            	width: 400px;
            	padding: 10px;
            	border: 3px solid rgba(0, 0, 0, 0.15);
            	border-radius: 10px;
            	box-shadow: 0 0 11px 0px #00000010;
            	/* That is, viewport minus top minus Live Docs */
            	max-height: calc(100vh - 5rem - 56px);
            	overflow: auto;
            	z-index: 40;
            	background-color: var(--main-bg-color);
            	transition: transform 300ms cubic-bezier(0.18, 0.89, 0.45, 1.12);

            			}
            			aside.plutoui-aside-wrapper > div {
            #				width: 300px;
            			}
            		}
            		</style>

            		<aside class="plutoui-aside-wrapper-$(uuid)">
            		<div>
            		$(text)
            		</div>
            		</aside>

            		"""
        )
    end
    floataside(stuff; kwargs...) = floataside(md"""$(stuff)"""; kwargs...)
end;

# ╔═╡ ad804cc3-93cd-421b-8314-a7d11e051e3d
floataside(
    md"""
    Select voltage in plots:

    $(@bind pb_v PlutoUI.Slider(range(0,vmax,length=201),default=0.05,show_value=true))
    """, top = 400
)

# ╔═╡ 811fca9b-bac0-4003-a0f9-bfefcbfbfa30
plotcdl(pbo_volts, pbo_caps, pb_v)

# ╔═╡ 65a25bd3-c420-4da3-b808-0e1888c6b4ba
plotsols(pbo_sols, pb_v, pbo_molfrac)

# ╔═╡ eb7d0765-d5e7-4ef9-916d-764c5aca9822
plotcdl(pbi_volts, pbi_caps, pb_v)

# ╔═╡ 5d15d24d-6317-4ad1-a53e-5af5c5bcf28a
plotsols(pbi_sols, pb_v, pbi_molfrac)

# ╔═╡ 4652157f-3288-45b6-877e-17ed8dc7d18e
function slidemodeswitch()
    uuid = uuid1()
    return html"""
    	<script id=$(uuid)>

        const right = document.querySelector('button.changeslide.next')
        const left = document.querySelector('button.changeslide.prev')

        let fullScreen = false

        const func = (e) => {
            if (e.key == "F10") {
                e.preventDefault()
                window.present()
                if (fullScreen) {
                    document.exitFullscreen().then(() => fullScreen = false)
                } else {

                    document.documentElement.requestFullscreen().then(() => fullScreen = true)
                }
            }
            if (document.body.classList.contains('presentation')) {
             
                if (e.target.tagName == "TEXTAREA") return
                if (e.key == "PageUp") {
                    e.preventDefault()
                    left.click()
                    return
             }

                if (e.key == "PageDown") {
                    e.preventDefault()
                    right.click()
                    return
                }
                if (e.key == "Escape") {
                    window.present()
                    fullScreen = false
                document.exitFullscreen().catch(() => {return})
                }
            }
        }

        document.addEventListener('keydown',func)

        invalidation.then(() => {document.removeEventListener('keydown',func)})
    </script>
    """
end;


# ╔═╡ Cell order:
# ╟─66dc2136-a601-4291-9aa5-81616a4efdd7
# ╟─c8ac41dd-df4d-47f6-909e-4e18e65ca994
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╟─b181a82b-46e7-4dae-b0b1-6886dd40886a
# ╟─5a5e3586-2009-4381-ae9d-a0bd34539b65
# ╠═b6338637-1f3c-4410-8e82-3afdaf603656
# ╟─371df98d-9f55-4244-8243-e6b597cc8d84
# ╠═10be79ba-9ab9-4e6f-8568-2de93ea77964
# ╟─834c81c8-9035-449f-b331-73bbee87756c
# ╠═7f209979-0776-40ab-9b2b-3b2d0145fa2c
# ╟─80a42b07-d6b1-41d6-bf1d-66a01eeb5019
# ╠═fb653542-7084-4e6f-87b9-9ab29e89ba04
# ╟─c57e82bb-79d8-4553-b31d-d1448c38c649
# ╟─e2179a6c-cc1c-4850-89d4-d3f87fd5e6ee
# ╠═7da889ce-9c6b-4abc-b19d-6311aacc32b1
# ╠═a2a8132a-d54d-48a7-aa02-ff83593f0e16
# ╠═db023dac-5ce0-4126-bfd0-9036ccec11f4
# ╟─7ee7be9f-038d-4080-af7c-01ae8f7d2682
# ╟─40282ce8-e032-4524-9a20-971826a6bb97
# ╟─1b13b664-9d5e-43a1-a7b2-2a1c6da4c947
# ╠═5befbcff-2e7c-48db-a327-27dbb6779d63
# ╟─91ecd9f8-5823-4ba0-9e39-a644f24424ef
# ╟─f8914d81-117d-4d9e-8f89-27a3045c8e7f
# ╠═2a7c5a4f-48b0-4995-a37b-70fd3d5bf853
# ╠═ddbb2a10-f22e-4474-b02d-e7b514852e58
# ╟─a4a4c6a4-ea90-4d23-b13a-2020790b2889
# ╠═286174fc-9a64-426a-9d3c-a96b39f7ea0c
# ╠═39996246-cb69-4c59-802a-d00804d84b0b
# ╠═95e49176-6c34-4c89-b156-5972080f82bb
# ╠═91a661ee-9c0c-496e-b6fb-4bb692cb7255
# ╠═0c77e9c2-1242-4694-bacf-3b6c9856997a
# ╠═163c5395-726f-4874-9304-ab9a46cc8bf9
# ╠═811fca9b-bac0-4003-a0f9-bfefcbfbfa30
# ╠═65a25bd3-c420-4da3-b808-0e1888c6b4ba
# ╟─db5b6820-5a53-465c-b380-66756fd722a6
# ╟─11bbb453-0a80-4ec4-b5a9-c044ebce8146
# ╠═c2500d41-53e8-4ac7-8e52-5107cf0335b2
# ╠═cb8b5a32-c414-4768-b250-2398cc432881
# ╠═77e92a09-66e8-4502-8115-c892eef1f597
# ╠═0b628db9-1719-4e1a-b87b-92d3e2808e5c
# ╠═db911ce1-c0f1-4261-aee7-54156ed03248
# ╟─eb7d0765-d5e7-4ef9-916d-764c5aca9822
# ╟─5d15d24d-6317-4ad1-a53e-5af5c5bcf28a
# ╟─9b5f389f-b105-4610-bab7-f79305fedc31
# ╟─ad804cc3-93cd-421b-8314-a7d11e051e3d
# ╟─e373824d-6e0d-44c8-8c8a-192cccf10298
# ╠═a290a533-9639-465f-b6ef-66c28d6136e8
# ╟─5beb3a0d-e57a-4aea-b7a0-59b8ce9ff5ce
# ╟─f9b4d4dc-7def-409f-b40a-f4eba1163741
# ╟─7a93e9a8-8a2d-4b11-84ef-691706c0eb0f
# ╟─86e1fffe-cf70-4d93-8131-153f0e98c94b
# ╟─4652157f-3288-45b6-877e-17ed8dc7d18e
