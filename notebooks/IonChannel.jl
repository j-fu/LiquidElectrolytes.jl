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
	bfacemask!(pnpgrid,[0,L/2], [W, L/2], 5 )
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

    boundary_dirichlet!(y, u, bnode, species = ip, region = 1, value = 0)
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
	    ip=4
	    iϕ=3
	    nrow=2
	    c1=sol[1, :] 
	    c2=sol[2, :]
	    grid1=grid
		if size(sol,1)>4
		    grid1=subgrid(grid,[1])
		    c1=view(sol[1, :],grid1) 
		    c2=view(sol[2, :],grid1) 
	
			grid2=subgrid(grid,[2])
		    c3=view(sol[3, :],grid2) 
		    c4=view(sol[4, :],grid2) 

			
			iϕ=5
			ip=6
			nrow=3
		end
		climits=(0,5)
        vis = GridVisualizer(; layout = (nrow, 2), size = (650, 600))
        scalarplot!(vis[1, 1], grid1, c1/ cbulk; limits=climits)
        scalarplot!(vis[1, 2], grid1, c2 / cbulk; limits=climits)
	    if size(sol,1)>4
	        scalarplot!(vis[2, 1], grid2, c3 / cbulk; limits=climits)
	        scalarplot!(vis[2, 2], grid2, c4 / cbulk; limits=climits)
		end
        scalarplot!(vis[nrow, 1], grid, sol[iϕ, :], colormap = :bwr, limits = (-1, 1))
        scalarplot!(vis[nrow, 2], grid, sol[ip, :])
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

# ╔═╡ 79a297c3-931e-4f4d-b0b0-27cde5f69dee
begin
	dgldata1=PNPData()
	dgldata2=PNPData()
	dgldata2.D*=0.1
end

# ╔═╡ f087896b-652f-48bc-ad78-9d66fa339f89
struct MyCellData{TE} <: AbstractCellData
	electrolytes::TE
end

# ╔═╡ faff62ea-e2df-4fd1-8872-fb0e98bf431f
LiquidElectrolytes.electrolytes(data::MyCellData)= data.electrolytes

# ╔═╡ e892d2f7-1f58-4d94-9719-bec80588ec97
mycelldata=MyCellData([dgldata1, dgldata2])

# ╔═╡ 19e2ff08-5347-41e1-bef8-808e415eb3ed
md"""
## Two electrolytes with interface reaction
"""

# ╔═╡ 57d2d467-6e2e-4501-be9e-7b2e643c8c16
Base.@kwdef struct MyCellData2{TE} <: AbstractCellData
	electrolytes::TE
	k::Float64=1.0e1
end

# ╔═╡ 771bf8d8-0e68-4ea5-b926-9fa945d239ba
LiquidElectrolytes.electrolytes(data::MyCellData2)= data.electrolytes

# ╔═╡ 20d53525-cbb3-44a4-bce2-d18fb395c66a
function pnpbcond2(y, u, bnode, data)
    (; iϕ, ip, cspecies) = electrolytes(data)[1]
	λ=bnode.embedparam
    boundary_neumann!(y, u, bnode, species = iϕ, region = 2, value = σ*λ)
    boundary_dirichlet!(y, u, bnode, species = iϕ, region = 1, value = -Δϕ / 2)
    boundary_dirichlet!(y, u, bnode, species = iϕ, region = 3, value = Δϕ / 2)
    boundary_dirichlet!(y, u, bnode, species = ip, region = 1, value = 0)
    boundary_dirichlet!(y, u, bnode, species = ip, region = 3, value = 0)
    for i in cspecies
        boundary_dirichlet!(y, u, bnode, species = i, region = 1, value = cbulk)
        boundary_dirichlet!(y, u, bnode, species = i, region = 3, value = cbulk)
    end
    return
end

# ╔═╡ d466c44b-48e4-4b17-bcfc-44f01cca90c5
sys2=PNPSystem(pnpgrid, celldata=mycelldata, bcondition=pnpbcond2)

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

# ╔═╡ 2d7242e7-4814-4192-af57-8ec3a2df18d9
edata1=ElectrolyteData(iϕ=5, ip=6, cspecies=[1,2], nc=4)

# ╔═╡ ec12834e-dad9-48fb-b6e7-7b99d2a75d49
edata2=ElectrolyteData(iϕ=5, ip=6, cspecies=[3,4], nc=4)

# ╔═╡ 51257300-35c3-4fc4-bbf0-0bae0ab2fc37
mycelldata2=MyCellData2(electrolytes=[edata1, edata2])

# ╔═╡ 65e04003-6e2b-4cdc-b09a-ad41195f21d6
function pnpbcond3(y, u, bnode, data)
	elytes=electrolytes(data)
	k=data.k
    (; iϕ, ip) = elytes[1]
	if bnode.region==5
		r13= k*(u[1] - u[3])
		y[1]+=r13
		y[3]-=r13
		r24= k*(u[2] - u[4])
		y[2]+=r24
		y[4]-=r24
	end

	
	λ=bnode.embedparam
    boundary_neumann!(y, u, bnode, species = iϕ, region = 2, value = σ*λ)
    boundary_dirichlet!(y, u, bnode, species = iϕ, region = 1, value = -λ*Δϕ / 2)
    boundary_dirichlet!(y, u, bnode, species = iϕ, region = 3, value = λ*Δϕ / 2)

	boundary_dirichlet!(y, u, bnode, species = ip, region = 1, value = 0)
    boundary_dirichlet!(y, u, bnode, species = ip, region = 3, value = 0)
    for i in elytes[1].cspecies
        boundary_dirichlet!(y, u, bnode, species = i, region = 3, value = cbulk)
    end
	for i in elytes[2].cspecies
        boundary_dirichlet!(y, u, bnode, species = i, region = 1, value = cbulk)
    end
    return
end

# ╔═╡ 85029a95-f297-4c53-b3b1-13b168d1b59a


# ╔═╡ c3acdde2-7f4c-415f-b8d7-d6f2b27edf92
sys3=PNPSystem(pnpgrid, celldata=mycelldata2, bcondition=pnpbcond3,
			  unknown_storage=:dense)

# ╔═╡ 60ce3f69-4d31-448b-9623-741947ea40ce
begin
		uini=unknowns(sys3.vfvmsys, inival=0)
	    elytes=electrolytes(mycelldata2)
		uini[1,:].= elytes[1].c_bulk[1]
		uini[2,:].= elytes[1].c_bulk[2]
		uini[3,:].= elytes[2].c_bulk[3]
		uini[4,:].= elytes[2].c_bulk[4]
end

# ╔═╡ 0bb425f7-f135-489c-81bc-9128fb87cd67
esol3=solve(sys3.vfvmsys; inival=uini, embed=[0,1], Δp=2.0e-1, Δp_min=1.0e-2 , 
			Δu_opt=1.0e4, verbose="ne", 
			damp_initial=1, force_first_step=true, 
			max_round=3, tol_round=1.0e-7, maxiters=100)

# ╔═╡ dca04b49-3d52-43b4-ad30-c90bcc26f74e
#=╠═╡
if isdefined(Main,:PlutoRunner)
	plotsol(pnpgrid, esol3[end])
end
  ╠═╡ =#

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
# ╠═f087896b-652f-48bc-ad78-9d66fa339f89
# ╠═faff62ea-e2df-4fd1-8872-fb0e98bf431f
# ╠═e892d2f7-1f58-4d94-9719-bec80588ec97
# ╠═d466c44b-48e4-4b17-bcfc-44f01cca90c5
# ╠═fff71cb0-c109-4686-8411-74528a98714c
# ╠═d5f33923-57f7-4e08-88fb-78c84c147a98
# ╠═2ccdf56a-6906-49af-81a3-00f32ba4e40f
# ╟─19e2ff08-5347-41e1-bef8-808e415eb3ed
# ╠═57d2d467-6e2e-4501-be9e-7b2e643c8c16
# ╠═771bf8d8-0e68-4ea5-b926-9fa945d239ba
# ╠═2d7242e7-4814-4192-af57-8ec3a2df18d9
# ╠═ec12834e-dad9-48fb-b6e7-7b99d2a75d49
# ╠═51257300-35c3-4fc4-bbf0-0bae0ab2fc37
# ╠═65e04003-6e2b-4cdc-b09a-ad41195f21d6
# ╠═85029a95-f297-4c53-b3b1-13b168d1b59a
# ╠═c3acdde2-7f4c-415f-b8d7-d6f2b27edf92
# ╠═60ce3f69-4d31-448b-9623-741947ea40ce
# ╠═0bb425f7-f135-489c-81bc-9128fb87cd67
# ╠═dca04b49-3d52-43b4-ad30-c90bcc26f74e
# ╟─8af12f1c-d35b-4cc9-8185-1bb5adbb69e8
# ╟─baa3e08e-5d64-4c8f-9f6d-5fdb40e97bc5
# ╟─784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
# ╟─afe4745f-f9f1-4e23-8735-cbec6fb79c41
