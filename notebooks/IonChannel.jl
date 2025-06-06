### A Pluto.jl notebook ###
# v0.20.9

using Markdown
using InteractiveUtils

# ╔═╡ a70cef7d-2a2f-4155-bdf3-fec9df94c63f
begin
    using LiquidElectrolytes
    using Printf
    using ExtendableGrids
    using LinearAlgebra
    using PlutoUI
    using VoronoiFVM
    using LessUnitful
    using Test

end

# ╔═╡ 4a215269-da12-431b-9e56-beec46112998
# ╠═╡ skip_as_script = true
#=╠═╡
begin
   using Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "docs")) # hide
    using Revise
	    using GridVisualize

   	import CairoMakie	
	    using CairoMakie: lines!
  	 	default_plotter!(CairoMakie)
 		CairoMakie.activate!(type="png")
    TableOfContents()
end
  ╠═╡ =#

# ╔═╡ baa3e08e-5d64-4c8f-9f6d-5fdb40e97bc5
begin
    using HypertextLiteral: @htl_str, @htl
    using UUIDs: uuid1
end

# ╔═╡ 10ecf504-a07c-4bbd-a340-fbba905738f4
md"""
## General  setup
"""

# ╔═╡ 8db5d72c-34a9-411e-ade5-5f48b99f9231
begin
    @phconstants e N_A
    const F = e * N_A
    @unitfactors mol m dm cm nm μA s mV V mPa

end

# ╔═╡ 8af5f374-d379-44a6-b567-ca34f6a2149e
begin
    const L = 20 * nm
    const W = 5 * nm
    const nref = 0
    const cyl = false
    nr = 8 * 2^nref + 1
    nl = Int((nr - 1) * L / W) + 1
    Y = range(0, L, length = nl)
    pnpX = geomspace(0, W, 2 * W / nr, 0.25 * W / nr)
    pnpgrid = simplexgrid(pnpX, Y)
	cellmask!(pnpgrid, [0,0], [W,L/2], 2)
end

# ╔═╡ f49f1197-257b-4fdb-ad80-5236c3aa1ee7
#=╠═╡
let
    if isdefined(Main,:PlutoRunner)
        vis = GridVisualizer(layout = (1, 2), size = (650, 300), linewidth = 0.1)
        gridplot!(vis[1, 2], pnpgrid, title = "pnpgrid", aspect = 0.5)
        reveal(vis)
    end
end
  ╠═╡ =#

# ╔═╡ 154c142a-c2b2-4229-80ec-9d7cf924c04c
begin
    const cbulk = 1.0 * mol / dm^3
    const Δϕ = 0.5 * V
    const σ = -15* μA * s / cm^2
end

# ╔═╡ 280f5233-621e-4d9d-95f2-8c97a2c343fc
default_pnpdata = ElectrolyteData();

# ╔═╡ eb84642a-a26d-4b26-8530-35a499a179a8
function PNPData(; molar_volume = default_pnpdata.v0, κ = 5)
    pnpdata = ElectrolyteData()
    pnpdata.κ .= κ
    pnpdata.D .= 1.0e-9 * m^2 / s
    pnpdata.v .= molar_volume
    update_derived!(pnpdata)
    return pnpdata
end

# ╔═╡ 93661e00-4533-4e4e-a9f4-ffa31feabd3b
md"""
## One electrolyte case
"""

# ╔═╡ 21d159e0-1639-451e-9054-2729a5942b5f
function pnpbcond(y, u, bnode, data)
    (; iϕ, ip, cspecies) = data
	λ=bnode.embedparam
    boundary_neumann!(y, u, bnode, species = iϕ, region = 2, value = σ*λ)
    boundary_dirichlet!(y, u, bnode, species = iϕ, region = 1, value = -Δϕ / 2)
    boundary_dirichlet!(y, u, bnode, species = iϕ, region = 3, value = Δϕ / 2)
    boundary_dirichlet!(y, u, bnode, species = ip, region = 3, value = 0)
    for i in cspecies
        boundary_dirichlet!(y, u, bnode, species = i, region = 1, value = cbulk)
        boundary_dirichlet!(y, u, bnode, species = i, region = 3, value = cbulk)
    end
    return
end

# ╔═╡ f78d5edc-9e67-4d93-87f2-a08a9fc3b845
dgldata = PNPData()

# ╔═╡ 259fa06a-fad4-4551-9e30-55eafdc98584
sys=PNPSystem(pnpgrid, celldata=PNPData(), bcondition=pnpbcond)

# ╔═╡ 837be65c-f2d8-4800-96e7-ee27d448bf24
inival=unknowns(sys)

# ╔═╡ 2b120621-9bf6-495b-8f27-9185a8c54c7c
esol=solve(sys.vfvmsys; inival, embed=[0,1], Δp=0.1, Δu_opt=1.0e3, verbose="ne")

# ╔═╡ bd329c64-a346-4692-a79c-b8b04106b785
dglpnpsol=esol[end]

# ╔═╡ 0e558929-d8b2-4c22-9596-322f0e6f335d
#=╠═╡
function plotsol(grid,sol)
	
        vis = GridVisualizer(; layout = (2, 2), size = (650, 600))
        scalarplot!(vis[1, 1], grid, sol[1, :] / cbulk)
        scalarplot!(vis[1, 2], grid, sol[2, :] / cbulk)
        scalarplot!(vis[2, 1], grid, sol[3, :], colormap = :bwr, limits = (-1, 1))
        scalarplot!(vis[2, 2], grid, sol[4, :])
        reveal(vis)
    end

  ╠═╡ =#

# ╔═╡ a91c5744-e744-4c25-a748-30d7155e9b63
#=╠═╡
if isdefined(Main,:PlutoRunner)
	plotsol(pnpgrid, esol[end])
end
  ╠═╡ =#

# ╔═╡ e0834845-b18f-4fbb-8780-6d8e09b31408
md"""
## Two electrolyte case
"""

# ╔═╡ 20d53525-cbb3-44a4-bce2-d18fb395c66a
function pnpbcond2(y, u, bnode, data)
    (; iϕ, ip, cspecies) = data[1]
	λ=bnode.embedparam
    boundary_neumann!(y, u, bnode, species = iϕ, region = 2, value = σ*λ)
    boundary_dirichlet!(y, u, bnode, species = iϕ, region = 1, value = -Δϕ / 2)
    boundary_dirichlet!(y, u, bnode, species = iϕ, region = 3, value = Δϕ / 2)
    boundary_dirichlet!(y, u, bnode, species = ip, region = 3, value = 0)
    for i in cspecies
        boundary_dirichlet!(y, u, bnode, species = i, region = 1, value = cbulk)
        boundary_dirichlet!(y, u, bnode, species = i, region = 3, value = cbulk)
    end
    return
end

# ╔═╡ 79a297c3-931e-4f4d-b0b0-27cde5f69dee
begin
	dgldata1=PNPData()
	dgldata2=PNPData()
	dgldata2.D*=0.1
end

# ╔═╡ d466c44b-48e4-4b17-bcfc-44f01cca90c5
sys2=PNPSystem(pnpgrid, celldata=[dgldata1, dgldata2], bcondition=pnpbcond2)

# ╔═╡ fff71cb0-c109-4686-8411-74528a98714c
inival2=unknowns(sys2)

# ╔═╡ d5f33923-57f7-4e08-88fb-78c84c147a98
esol2=solve(sys2.vfvmsys; inival2, embed=[0,1], Δp=0.1, Δu_opt=1.0e3, verbose="ne")

# ╔═╡ 2ccdf56a-6906-49af-81a3-00f32ba4e40f
#=╠═╡
if isdefined(Main,:PlutoRunner)
	plotsol(pnpgrid, esol2[end])
end
  ╠═╡ =#

# ╔═╡ 19e2ff08-5347-41e1-bef8-808e415eb3ed
md"""
## Two electrolytes with interface reaction
"""

# ╔═╡ 8af12f1c-d35b-4cc9-8185-1bb5adbb69e8
html"""<hr>"""

# ╔═╡ 784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
html"""<style>.dont-panic{ display: none }</style>"""

# ╔═╡ afe4745f-f9f1-4e23-8735-cbec6fb79c41
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


# ╔═╡ Cell order:
# ╠═4a215269-da12-431b-9e56-beec46112998
# ╠═a70cef7d-2a2f-4155-bdf3-fec9df94c63f
# ╟─10ecf504-a07c-4bbd-a340-fbba905738f4
# ╠═8db5d72c-34a9-411e-ade5-5f48b99f9231
# ╠═8af5f374-d379-44a6-b567-ca34f6a2149e
# ╠═f49f1197-257b-4fdb-ad80-5236c3aa1ee7
# ╠═154c142a-c2b2-4229-80ec-9d7cf924c04c
# ╠═280f5233-621e-4d9d-95f2-8c97a2c343fc
# ╠═eb84642a-a26d-4b26-8530-35a499a179a8
# ╟─93661e00-4533-4e4e-a9f4-ffa31feabd3b
# ╠═21d159e0-1639-451e-9054-2729a5942b5f
# ╠═f78d5edc-9e67-4d93-87f2-a08a9fc3b845
# ╠═259fa06a-fad4-4551-9e30-55eafdc98584
# ╠═837be65c-f2d8-4800-96e7-ee27d448bf24
# ╠═2b120621-9bf6-495b-8f27-9185a8c54c7c
# ╠═bd329c64-a346-4692-a79c-b8b04106b785
# ╠═0e558929-d8b2-4c22-9596-322f0e6f335d
# ╠═a91c5744-e744-4c25-a748-30d7155e9b63
# ╟─e0834845-b18f-4fbb-8780-6d8e09b31408
# ╠═20d53525-cbb3-44a4-bce2-d18fb395c66a
# ╠═79a297c3-931e-4f4d-b0b0-27cde5f69dee
# ╠═d466c44b-48e4-4b17-bcfc-44f01cca90c5
# ╠═fff71cb0-c109-4686-8411-74528a98714c
# ╠═d5f33923-57f7-4e08-88fb-78c84c147a98
# ╠═2ccdf56a-6906-49af-81a3-00f32ba4e40f
# ╟─19e2ff08-5347-41e1-bef8-808e415eb3ed
# ╟─8af12f1c-d35b-4cc9-8185-1bb5adbb69e8
# ╟─baa3e08e-5d64-4c8f-9f6d-5fdb40e97bc5
# ╟─784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
# ╟─afe4745f-f9f1-4e23-8735-cbec6fb79c41
