### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 4fc7fda6-423b-48ea-8f86-6718a9050ee0
begin
    import Pkg # hide
    Pkg.activate(joinpath(@__DIR__, "..", "docs")) # hide
    using Revise # hide
    using PlutoUI, HypertextLiteral, UUIDs, Markdown
    using Test
    using LinearAlgebra
    using LessUnitful
    using VoronoiFVM
    using LiquidElectrolytes
    using ExtendableGrids
    using LessUnitful
    using GridVisualize
    if isdefined(Main, :PlutoRunner)
        using Colors
        using CairoMakie
        default_plotter!(CairoMakie)
        CairoMakie.activate!(type = "png")
    end
end;

# ╔═╡ 9391dc5a-9249-4719-9347-98c95c80499d
md"""
# DLCap.jl

Comparison of equilibrium results for different flux and equilibrium implementations

[Pluto source](https://raw.githubusercontent.com/j-fu/LiquidElectrolytes.jl/main/notebooks/DLCap.jl)
"""

# ╔═╡ e04f6d4c-144c-4de2-83b4-11190ad0b07d
TableOfContents(aside = false)

# ╔═╡ f5c18bb8-70d5-4cb4-88c8-4333d539d08e
pkgdir(LiquidElectrolytes)

# ╔═╡ 8c4a3420-191d-41e5-abbb-a50f2e16cd1c
begin
    const nm = ufac"nm"
    const V = ufac"V"
    const mol = ufac"mol"
    const dm = ufac"dm"
end;

# ╔═╡ a5193642-e9fc-4b80-a9eb-84f5e1795008
begin
    const molarity = 0.01
    const voltages = -1:0.01:1
    const nref = 1
    const solvation = 10
    celldata = ElectrolyteData(;
        nc = 2,
        Γ_we = 1,
        Γ_bulk = 2
    )

    celldata.c_bulk .= molarity * mol / dm^3
    celldata.κ .= solvation
    celldata
end


# ╔═╡ 4a07ae69-5db0-4868-a0b1-ee4b22e95eb6
grid = let
    hmin = 1.0e-1 * nm * 2.0^(-nref)
    hmax = 1.0 * nm * 2.0^(-nref)
    L = 20.0 * nm
    X = geomspace(0, L, hmin, hmax)
    simplexgrid(X)
end

# ╔═╡ 229ae84b-fcff-4767-b03b-c17b465de57d
gridplot(grid, size = (600, 200))

# ╔═╡ 375af0c8-b008-4f8a-9cfb-ba10f2c2787b
## Define boundary conditions
function bcondition(f, u, bnode, data)
    (; Γ_we, Γ_bulk, ϕ_we) = data
    iϕ, ip = 1, 2
    ## Dirichlet ϕ=ϕ_we at Γ_we
    boundary_dirichlet!(f, u, bnode, species = iϕ, region = Γ_we, value = ϕ_we)
    boundary_dirichlet!(f, u, bnode, species = iϕ, region = Γ_bulk, value = data.ϕ_bulk)
    return boundary_dirichlet!(f, u, bnode, species = ip, region = Γ_bulk, value = data.p_bulk)
end


# ╔═╡ 0ba45da2-b5e7-4eb9-ab37-2bcfbe1307f4
pbsys = PBSystem(grid; celldata = deepcopy(celldata), bcondition)

# ╔═╡ 2ad3f1ad-3bc9-4ced-a9d8-11c770f9710f
result = dlcapsweep(
    pbsys;
    inival = unknowns(pbsys),
    iϕ = 1, voltages,
    molarity = molarity * ufac"mol/dm^3",
    damp_initial = 0.1
)

# ╔═╡ 6b0adcfd-a279-4b58-b029-db07fca60b0f
begin
    ecell = create_equilibrium_system(grid, EquilibriumData(celldata))
    evolts, ecaps = dlcapsweep_equi(
        ecell;
        vmax = voltages[end],
        molarity,
        δV = 1.0e-4,
        nsteps = 101
    )
end;

# ╔═╡ f57cb1f9-3505-45e1-9524-ba37af47637b
function pnpsweep(scheme)

    function pnp_bcondition(f, u, bnode, data::ElectrolyteData)
        (; iϕ, Γ_we, ϕ_we) = data
        boundary_dirichlet!(f, u, bnode, species = iϕ, region = Γ_we, value = ϕ_we)
        return bulkbcondition(f, u, bnode, data)
    end
    xcelldata = deepcopy(celldata)
    xcelldata.scheme = scheme
    μcell = PNPSystem(grid; bcondition = pnp_bcondition, celldata = xcelldata)
    return μresult = dlcapsweep(
        μcell;
        voltages,
        molarity = molarity * ufac"mol/dm^3",
        δ = 1.0e-4, damp_initial = 0.1
    )
end;

# ╔═╡ 8e182773-b012-446f-b676-39d96e2ea424
μresult = pnpsweep(:μex)

# ╔═╡ 284acfe8-5c67-4704-b055-78f1acba4d76
aresult = pnpsweep(:act)

# ╔═╡ aaa16bc5-a00d-4172-a91c-04f72a653938
cresult = pnpsweep(:cent)

# ╔═╡ 0ce644b5-4375-4157-8075-c1beac3cb789
myround(x) = round(x, sigdigits = 5)

# ╔═╡ 38a6127d-b806-4056-bd60-aac5bb717a4e
@test isapprox(result.dlcaps, ecaps, rtol = 5.0e-4)

# ╔═╡ 43d18b49-7107-49b6-8073-ffb6dc39dcea
@test isapprox(result.dlcaps, μresult.dlcaps, rtol = 1.0e-11)

# ╔═╡ d59914df-9b6e-4448-84af-8a24750ee0a0
@test isapprox(result.dlcaps, aresult.dlcaps, rtol = 1.0e-11)

# ╔═╡ fa6dfae5-4b92-4e51-b55f-d6b3d1180ee7
@test isapprox(result.dlcaps, cresult.dlcaps, rtol = 1.0e-11)

# ╔═╡ ab9fa489-e718-4b2c-8689-cb627fd30593
typeof(md""" abc""") |> supertype

# ╔═╡ 7a93e9a8-8a2d-4b11-84ef-691706c0eb0f
begin
    hrule() = html"""<hr>"""
    highlight(mdstring, color) = htl"""<blockquote style="padding: 10px; background-color: $(color);">$(mdstring)</blockquote>"""

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

# ╔═╡ 5beb3a0d-e57a-4aea-b7a0-59b8ce9ff5ce
hrule()

# ╔═╡ b3b4a9bb-022e-4cd4-8547-b1b3d6176af2
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
end

# ╔═╡ 222b34e1-8d3c-4143-8826-7b522c19e23c
floataside(
    let
        vis = GridVisualizer(size = (350, 200), legend = :rt)
        scalarplot!(vis, result.voltages, result.dlcaps, color = :red, label = "pb")

        scalarplot!(vis, evolts, ecaps, color = :green, label = "equi", clear = false)
        scalarplot!(vis, μresult.voltages, μresult.dlcaps, color = :blue, label = "μpnp", clear = false)
        scalarplot!(vis, aresult.voltages, aresult.dlcaps, color = :magenta, label = "apnp", clear = false)
        scalarplot!(vis, cresult.voltages, cresult.dlcaps, color = :cyan, label = "cpnp", clear = false)

        reveal(vis)
    end, top = 200
)

# ╔═╡ 948b0ebb-87ee-4c1e-b854-cb13db390934
let
    pb0 = findfirst(x -> x ≈ 0, result.voltages)
    equi0 = findfirst(x -> x ≈ 0, evolts)
    pnp0 = findfirst(x -> x ≈ 0, μresult.voltages)

    floataside(
        md"""
        DLcap0: $(dlcap0(celldata)|>myround), $(dlcap0(EquilibriumData(celldata))|>myround)

        DLCapMin(pb,equi,pnp): $(result.dlcaps[pb0]|>myround,ecaps[equi0]|>myround,μresult.dlcaps[pnp0]|>myround))
        	
        	""",

        top = 450

    )
end


# ╔═╡ Cell order:
# ╠═4fc7fda6-423b-48ea-8f86-6718a9050ee0
# ╟─9391dc5a-9249-4719-9347-98c95c80499d
# ╟─e04f6d4c-144c-4de2-83b4-11190ad0b07d
# ╠═f5c18bb8-70d5-4cb4-88c8-4333d539d08e
# ╠═8c4a3420-191d-41e5-abbb-a50f2e16cd1c
# ╠═4a07ae69-5db0-4868-a0b1-ee4b22e95eb6
# ╠═229ae84b-fcff-4767-b03b-c17b465de57d
# ╠═a5193642-e9fc-4b80-a9eb-84f5e1795008
# ╠═375af0c8-b008-4f8a-9cfb-ba10f2c2787b
# ╠═0ba45da2-b5e7-4eb9-ab37-2bcfbe1307f4
# ╠═2ad3f1ad-3bc9-4ced-a9d8-11c770f9710f
# ╠═6b0adcfd-a279-4b58-b029-db07fca60b0f
# ╠═f57cb1f9-3505-45e1-9524-ba37af47637b
# ╠═8e182773-b012-446f-b676-39d96e2ea424
# ╠═284acfe8-5c67-4704-b055-78f1acba4d76
# ╠═aaa16bc5-a00d-4172-a91c-04f72a653938
# ╠═0ce644b5-4375-4157-8075-c1beac3cb789
# ╠═222b34e1-8d3c-4143-8826-7b522c19e23c
# ╠═38a6127d-b806-4056-bd60-aac5bb717a4e
# ╠═43d18b49-7107-49b6-8073-ffb6dc39dcea
# ╠═d59914df-9b6e-4448-84af-8a24750ee0a0
# ╠═fa6dfae5-4b92-4e51-b55f-d6b3d1180ee7
# ╠═948b0ebb-87ee-4c1e-b854-cb13db390934
# ╟─5beb3a0d-e57a-4aea-b7a0-59b8ce9ff5ce
# ╟─ab9fa489-e718-4b2c-8689-cb627fd30593
# ╟─7a93e9a8-8a2d-4b11-84ef-691706c0eb0f
# ╟─b3b4a9bb-022e-4cd4-8547-b1b3d6176af2
