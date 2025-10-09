### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ a70cef7d-2a2f-4155-bdf3-fec9df94c63f
begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
    using Revise
	using PlutoUI
    using LessUnitful
    using LaTeXStrings
    using VoronoiFVM, ExtendableGrids, GridVisualize
    using LiquidElectrolytes
    using CairoMakie
    using Colors
    using DoubleFloats
    CairoMakie.activate!(; type = "png", visible = false)
    GridVisualize.default_plotter!(CairoMakie)
	TableOfContents()
end

# ╔═╡ 4219cd3c-8d65-4054-bdce-01410edee958
md"""
# Generate the figures in description.tex
"""

# ╔═╡ f674a3f5-87bf-4677-9081-4fd7b0cbb6e3
begin
    @phconstants N_A e R ε_0
    const F = N_A * e
    @unitfactors cm μF mol dm s mA A nm V M m
end

# ╔═╡ 7022555d-e74b-4900-8636-5bd3c4e26fc2
begin
    rlog = RLog()
    #LiquidElectrolytes.rlog(x::Number)=rlog(x)
end

# ╔═╡ b3f2b24f-2d3a-4095-94e7-1f858d7e1519
md"""
## Double layer capacitance
"""

# ╔═╡ 2e6524ce-86f0-47cb-a344-c014638b0336
function bcondition(f, u, bnode, data)
    (; iϕ, Γ_we, Γ_bulk, ϕ_we) = data

    ## Dirichlet ϕ=ϕ_we at Γ_we
    boundary_dirichlet!(f, u, bnode, species = iϕ, region = Γ_we, value = ϕ_we)

    ## Bulk condition at Γ_bulk
    return bulkbcondition(f, u, bnode, data, region = Γ_bulk)
end


# ╔═╡ b19bc908-1c88-4906-a9f7-3b240d9805b0
function dlcap_calc(;
        voltages = -1:0.01:1,
        molarities = [0.001, 0.01, 0.1, 1],
        κ = [10, 10, 10, 10],
        nref = 0
    )
    results = []
    cdl0 = []
    hmin = 1.0e-1 * nm * 2.0^(-nref)
    hmax = 1.0 * nm * 2.0^(-nref)
    L = 20.0 * nm
    X = geomspace(0, L, hmin, hmax)
    grid = simplexgrid(X)
    celldata = ElectrolyteData(;
        nc = 2,
        Γ_we = 1,
        Γ_bulk = 2,
        c_bulk = fill(0.01M, 2)
    )
    cell = PNPSystem(grid; bcondition, celldata)
    for imol in 1:length(molarities)
        celldata.κ .= κ[imol]
        celldata.c_bulk .= molarities[imol] * M
        result = dlcapsweep(
            cell;
            δ = 1.0e-6,
            voltages = collect(voltages) * V,
            verbose = "",
            store_solutions = true
        )
        push!(cdl0, dlcap0(celldata))

        push!(results, result)
    end
    return (; results, molarities, voltages, κ, cdl0, X)
end

# ╔═╡ e55f583a-a249-4cc1-b310-bfde0cab9d12
function dlcap_plot(res, output)
    vis = GridVisualizer(;
        resolution = (500, 400),
        legend = :rt,
        clear = true,
        xlabel = L"\phi/V",
        ylabel = L"C_{dl}/(μF/cm^2)",
    )
    (; results, molarities, voltages, κ, cdl0) = res

    for imol in 1:length(results)

        color = RGB(1 - imol / length(molarities), 0, imol / length(molarities))

        result = results[imol]


        scalarplot!(
            vis, result.voltages / V, result.dlcaps / (μF / cm^2);
            color,
            clear = false,
            label = "$(molarities[imol])M, κ=$(κ[imol])",
            markershape = :none
        )

        scalarplot!(
            vis, [0], [cdl0[imol]] / (μF / cm^2);
            clear = false,
            markershape = :circle,
            markersize = 10,
            label = ""
        )
    end
    p = reveal(vis)
    GridVisualize.save(joinpath(@__DIR__, output), vis)
    return p
end

# ╔═╡ 8f429eee-57c8-476f-8b79-c4a8035ba0b4
md"""
### Molarity variation
"""

# ╔═╡ a84ec6c7-6ee1-46e7-8a19-6a4c355fce6f
molvar_results = dlcap_calc(;
    nref = 1,
    molarities = [0.001, 0.01, 0.1, 1],
    κ = [10.0, 10, 10, 10]
)

# ╔═╡ f470297d-2649-424b-abe2-8418120a779a
dlcap_plot(molvar_results, "example_dlcap_molvar.png")

# ╔═╡ 4f6680d5-de48-4ce1-93c2-8b646423f677
md"""
### Nanoscale/Device scale
"""

# ╔═╡ 283c2ae6-41ac-40e1-9d7d-7d3da85b6494
let

    iresult = 4
    res = molvar_results.results[iresult]
    mity = molvar_results.molarities[iresult]
    ivolt = 160
    volt = res.voltages[ivolt]
    @show iresult, mity
    @show ivolt, volt
    sol = molvar_results.results[iresult].solutions[ivolt]
    X = molvar_results.X / nm
    fig = Figure(size = (750, 500))

	xt=[-1,0,2,3,4,6,8]
	xtnames=Any[s for s in string.(xt)]
	xtnames[4]=L"x_{DL}"
    axl = Axis(
        fig[1, 1], ylabel = L"ϕ/V", xlabel = "x/nm",
        width = 300, height = 150, title = "Nanoscale",
		xticks=(xt, xtnames  )
    )
    axr = Axis(
        fig[1, 1], yaxisposition = :right, ylabel = L"c/(mol/dm^3)",
        width = 300, height = 150,
			xticks=(xt, xtnames )
    )


    xlims!(axr, (-1, 8))
    xlims!(axl, (-1, 8))
    ylims!(axr, (-1, 6))
    ylims!(axl, (-0.1, 0.6))

    poly!(axl, Point2f[(-1, -1), (0, -1), (0, 6), (-1, 6)], color = RGBA(0, 0, 1, 0.25))
    poly!(axl, Point2f[(0, -1), (3, -1), (3, 6), (0, 6)], color = RGBA(1, 0, 0, 0.25))
    poly!(
        axl, Point2f[(3, -1), (8, -1), (8, 6), (3, 6)],
        color = RGBA(0.6, 0.6, 0.6, 0.25)
    )


    lines!(axl, [-1, 0], [sol[3, 1], sol[3, 1]], color = :green, linewidth = 3)
    data = [
        lines!(axl, X, sol[3, :], color = :green, linewidth = 3)
        lines!(axr, X, sol[1, :] / (mol / dm^3), color = :red, linewidth = 3)
        lines!(axr, X, sol[2, :] / (mol / dm^3), color = :blue, linewidth = 3)
    ]
    axislegend(axl, data, [L"ϕ", L"c^+", L"c^-"], position = :rt, backgroundcolor = :transparent)

    text!(axl, -0.5, 0.0, text = "metal", justification = :left, rotation = π / 2)
    text!(axl, 1.5, 0.0, text = "double layer", justification = :left, rotation = π / 2)
    text!(axl, 3.5, 0.12, text = "electroneutral bulk", justification = :left)

    ax2l = Axis(
        fig[2, 1], ylabel = L"ϕ/V", xlabel = "x/mm",
        width = 600, height = 150, title = "Experimental cell scale"
    )
    ax2r = Axis(
        fig[2, 1], yaxisposition = :right, ylabel = L"c/(mol/dm^3)",
        width = 600, height = 150
    )
    xlims!(ax2r, (-1, 8))
    xlims!(ax2l, (-1, 8))
    ylims!(ax2r, (-1, 6))
    ylims!(ax2l, (-0.1, 0.6))

    poly!(ax2l, Point2f[(-1, -1), (-0.05, -1), (-0.05, 6), (-1, 6)], color = RGBA(0, 0, 1, 0.25))
    poly!(ax2l, Point2f[(-0.05, -1), (0.05, -1), (0.05, 6), (-0.05, 6)], color = RGBA(1, 0, 0, 0.25))
    poly!(ax2l, Point2f[(0.05, -1), (10, -1), (10, 6), (0.05, 6)], color = RGBA(0.6, 0.6, 0.6, 0.25))


    lines!(ax2l, [-1, -0.05], [sol[3, 1], sol[3, 1]], color = :green, linewidth = 3)
    text!(ax2l, -0.5, 0.0, text = "metal", justification = :left, rotation = π / 2)
    text!(ax2l, 3.5, 0.2, text = "electroneutral bulk", justification = :left)

    data = [
        lines!(ax2l, [0.05, 8], [0, 0], color = :green, linewidth = 3)
        lines!(ax2r, [0.05, 8], [1, 1], color = :red, linewidth = 3)
        lines!(ax2r, [0.05, 8], [1, 1], color = :blue, linewidth = 3)
    ]
    axislegend(ax2l, data, [L"ϕ", L"c^+", L"c^-"], position = :rt, backgroundcolor = :transparent)
    CairoMakie.save(joinpath(@__DIR__, "bvscales.png"), fig)


    fig


end

# ╔═╡ 59404b67-8af9-4650-8df4-b06dda355e11
md"""
### Solvation number variation
"""

# ╔═╡ ad85052b-f5cd-451b-a978-7df924b1cd61
kappavar_results = dlcap_calc(;
    molarities = [0.1, 0.1, 0.1, 0.1],
    κ = [5, 10, 15, 20], nref = 1
)

# ╔═╡ 17ba9584-d7a2-4745-94f7-5bd651e67fa5
dlcap_plot(kappavar_results, "example_dlcap_kappavar.png")

# ╔═╡ aef407a4-b3d4-403b-99b0-6fc29883ead9
md"""
## BVCompare
"""

# ╔═╡ a1b85df4-954f-4c77-9d07-ab6e0f6d5745
function bvcompare_calc(;
        nref = 0, zR = 1, n = 1, verbose = "", bgmol = 2,
        rlog = Base.log,
        valuetype = Float64,
        kwargs...
    )
    voltages = (-1:0.025:1) * ufac"V"
    dlcap = false
    molR = 0.1
    molO = 0.1
    molS = bgmol
    iO = 1
    iR = 2
    iS = 3
    iX = 4
    zO = zR + n
    zS = -1
    zX = 1
    function halfcellbc(f, u, bnode, data)
        (; nc, Γ_we, Γ_bulk, ϕ_we, ip, iϕ, v, v0, RT) = data
        bulkbcondition(f, u, bnode, data; region = Γ_bulk)
        if bnode.region == Γ_we
            if !data.eneutral
                boundary_dirichlet!(
                    f, u, bnode;
                    species = iϕ, region = Γ_we, value = ϕ_we
                )
            end
            c0, barc = c0_barc(u, data)
            μR = chemical_potential(u[iR], barc, u[ip], v[iR] + κ * v0, data)
            μO = chemical_potential(u[iO], barc, u[ip], v[iO] + κ * v0, data)
            A = (μR - μO + Δg + data.eneutral * (zR - zO) * F * (u[iϕ] - ϕ_we)) / RT
            r = rrate(R0, β, A)
            f[iR] -= r
            f[iO] += r
        end
        return nothing
    end
    molX = molS - zR * molR - zO * molO
    z = [zO, zR, zS, zX]
    c_bulk = [molO, molR, molS, molX] * mol / dm^3
    @assert abs(c_bulk' * z < 1.0e-12)
    scheme = :μex
    κ = 5.0
    R0 = 1.0e-10 * mol / (cm^2 * s)
    Δg = 0.0
    β = 0.5
    hmin = 1.0e-2 * nm * 2.0^(-nref)
    hmax = 1.0 * nm * 2.0^(-nref)
    L = 20.0 * nm
    X = geomspace(0, L, hmin, hmax)
    grid = simplexgrid(X)
    pnpdata = ElectrolyteData(;
        nc = 4,
        z, c_bulk,
        κ = fill(κ, 4),
        Γ_we = 1,
        Γ_bulk = 2,
        rlog
    )


    function sweep(; eneutral = true, tunnel = false, bikerman = true)
        celldata = deepcopy(pnpdata)
        celldata.eneutral = eneutral
        vtype = valuetype
        if eneutral
            vtype = Float64
        end
        if !bikerman
            celldata.v .= 0.0
            celldata.κ .= 0.0
        end

        if tunnel
            pnpcell = PNPSystem(
                grid;
                bcondition = tunnelbc,
                reaction = tunnelreaction,
                celldata, valuetype = vtype, kwargs...
            )
        else
            pnpcell = PNPSystem(
                grid;
                bcondition = halfcellbc,
                celldata, valuetype = vtype, kwargs...
            )
        end
        result = ivsweep(
            pnpcell; voltages,
            store_solutions = true,
            verbose, Δu_opt = 1.0e10
        )
        currs = -zR * LiquidElectrolytes.currents(result, iR, electrode = :bulk) -
            zO * LiquidElectrolytes.currents(result, iO, electrode = :bulk)
        sol = LiquidElectrolytes.voltages_solutions(result)
        return result, currs, sol
    end

    nnpresult, nnpcurrs, nnpsol = sweep(eneutral = true, tunnel = false)
    pnpresult, pnpcurrs, pnpsol = sweep(eneutral = false, tunnel = false)
    return (; pnpresult, pnpcurrs, pnpsol, nnpresult, nnpcurrs, nnpsol, zR, n)
end


# ╔═╡ 21ab6d00-e979-4fce-bf6a-9cd4995d86f3
function bvcompare_plot(res, output = nothing; scaling = [1, 100, 1.0e5], legend = :lt)
    (; pnpresult, pnpcurrs, pnpsol, nnpresult, nnpcurrs, nnpsol, zR, n) = res
    tunnel = false

    function yplot(fg, imax; xlabel = "", legend = nothing)
        ax = Axis(fg; ylabel = L"I/(mA/cm^2)", xlabel)
        ylims!(ax, [-imax, imax])
        lines!(
            ax,
            pnpresult.voltages, -pnpcurrs / (mA / cm^2),
            color = :red, label = "PNP"
        )

        lines!(
            ax,
            nnpresult.voltages, -nnpcurrs / (mA / cm^2),
            color = :green, label = "NNP",
            linestyle = :dot,
            linewidth = 4
        )

        if tunnel
            lines!(
                p, tresult.voltages, -tcurrs / (mA / cm^2),
                color = :blue, label = "PNP-tunnel"
            )
        end
        if !isnothing(legend)
            axislegend(;
                position = legend, labelsize = 10,
                backgroundcolor = RGBA(1.0, 1.0, 1.0, 0.5)
            )
        end
        return ax
    end
    fig = Figure(size = (600, 400), fontsize = 12, linewidth = 2.5)
    title = tunnel ? L"z_R=%$(zR), n=%$(n), β=%$(btxt)/cm" : L"z_R=%$(zR), n=%$(n)"
    yplot(fig[1, :], scaling[1]; legend)
    yplot(fig[2, :], scaling[2])
    yplot(fig[3, :], scaling[3]; xlabel = L"U/V")
    Label(fig[0, :], title, tellwidth = false, fontsize = 16)

    if !isnothing(output)
        CairoMakie.save(joinpath(@__DIR__, output), fig)
    end
    return fig

end

# ╔═╡ 98de0d01-9d90-484a-849e-f544ad382aa0
md"""
### z_R=-1
"""

# ╔═╡ 1d6728f1-c88f-48ab-85b4-afa9372bf1a3
bvc_resultm1 = bvcompare_calc(; zR = -1, n = 1)

# ╔═╡ 41a5b4ce-b793-4dd0-9a82-221abc12c5c5
bvcompare_plot(bvc_resultm1, "bvc_zr=-1.png", scaling = [1.5, 50, 1.0e5], legend = :rt)

# ╔═╡ 008cff1e-93bf-4e3b-bde2-2cc8f9de311c
md"""
### z_R=0
"""

# ╔═╡ 45165906-d37a-4cec-8c98-2276f46054e3
bvc_result0 = bvcompare_calc(; zR = 0, n = 1)

# ╔═╡ 16539832-295d-49ef-a05c-be11d5d1b1a4
bvcompare_plot(bvc_result0, "bvc_zr=0.png", scaling = [1.5, 50, 1.0e5])

# ╔═╡ 082d2262-4af8-43a0-be91-39ab5004be27
md"""
### z_R=1
"""

# ╔═╡ c58a3109-8eaf-4c55-86d4-f49158405b2f
bvc_resultp1 = bvcompare_calc(;
    zR = 1, n = 1,
    valuetype = BigFloat,
    #							 rlog=RLog(Double64)
)

# ╔═╡ 7080aff5-3f11-4c5d-8cfe-2f7995e95450
bvcompare_plot(bvc_resultp1, "bvc_zr=1.png", scaling = [0.4, 150, 1.0e5], legend = :lt)

# ╔═╡ 955a53ac-5a2f-463f-bf10-ca50159861b7
bvc_resultx = bvcompare_calc(; zR = -1, n = 2, bgmol = 0.1)

# ╔═╡ 920a3fc6-91d4-42f1-ad49-a446703b6c46
bvcompare_plot(bvc_resultx, "bvc_x.png", scaling = [0.5, 300, 5.0e5])

# ╔═╡ 8af12f1c-d35b-4cc9-8185-1bb5adbb69e8
html"""<hr>"""

# ╔═╡ 784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
html"""<style>.dont-panic{ display: none }</style>"""

# ╔═╡ Cell order:
# ╟─4219cd3c-8d65-4054-bdce-01410edee958
# ╠═a70cef7d-2a2f-4155-bdf3-fec9df94c63f
# ╠═f674a3f5-87bf-4677-9081-4fd7b0cbb6e3
# ╠═7022555d-e74b-4900-8636-5bd3c4e26fc2
# ╟─b3f2b24f-2d3a-4095-94e7-1f858d7e1519
# ╠═2e6524ce-86f0-47cb-a344-c014638b0336
# ╠═b19bc908-1c88-4906-a9f7-3b240d9805b0
# ╠═e55f583a-a249-4cc1-b310-bfde0cab9d12
# ╟─8f429eee-57c8-476f-8b79-c4a8035ba0b4
# ╠═a84ec6c7-6ee1-46e7-8a19-6a4c355fce6f
# ╠═f470297d-2649-424b-abe2-8418120a779a
# ╟─4f6680d5-de48-4ce1-93c2-8b646423f677
# ╠═283c2ae6-41ac-40e1-9d7d-7d3da85b6494
# ╟─59404b67-8af9-4650-8df4-b06dda355e11
# ╠═ad85052b-f5cd-451b-a978-7df924b1cd61
# ╠═17ba9584-d7a2-4745-94f7-5bd651e67fa5
# ╟─aef407a4-b3d4-403b-99b0-6fc29883ead9
# ╠═a1b85df4-954f-4c77-9d07-ab6e0f6d5745
# ╠═21ab6d00-e979-4fce-bf6a-9cd4995d86f3
# ╟─98de0d01-9d90-484a-849e-f544ad382aa0
# ╠═1d6728f1-c88f-48ab-85b4-afa9372bf1a3
# ╠═41a5b4ce-b793-4dd0-9a82-221abc12c5c5
# ╟─008cff1e-93bf-4e3b-bde2-2cc8f9de311c
# ╠═45165906-d37a-4cec-8c98-2276f46054e3
# ╠═16539832-295d-49ef-a05c-be11d5d1b1a4
# ╟─082d2262-4af8-43a0-be91-39ab5004be27
# ╠═c58a3109-8eaf-4c55-86d4-f49158405b2f
# ╠═7080aff5-3f11-4c5d-8cfe-2f7995e95450
# ╠═955a53ac-5a2f-463f-bf10-ca50159861b7
# ╠═920a3fc6-91d4-42f1-ad49-a446703b6c46
# ╟─8af12f1c-d35b-4cc9-8185-1bb5adbb69e8
# ╟─784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
