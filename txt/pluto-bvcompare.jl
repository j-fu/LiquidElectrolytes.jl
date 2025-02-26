### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ a70cef7d-2a2f-4155-bdf3-fec9df94c63f
begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "docs"))
    using Revise
    using LessUnitful
    using LaTeXStrings
    using VoronoiFVM, ExtendableGrids, GridVisualize
    using LiquidElectrolytes
    using CairoMakie
    using Colors
    CairoMakie.activate!(; type = "png", visible = false)
    GridVisualize.default_plotter!(CairoMakie)
end

# ╔═╡ f674a3f5-87bf-4677-9081-4fd7b0cbb6e3
begin
    @phconstants N_A e R ε_0
    const F = N_A * e
    @unitfactors cm μF mol dm s mA A nm V M m
end

# ╔═╡ b3f2b24f-2d3a-4095-94e7-1f858d7e1519
md"""
## Dlcap
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
        voltages = -2:0.01:2,
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
        result = dlcapsweep(
            cell;
            δ = 1.0e-6,
            voltages = collect(voltages) * V,
            molarity = molarities[imol] * M,
            verbose = ""
        )
        push!(cdl0, dlcap0(celldata))

        push!(results, result)
    end
    return (; results, molarities, voltages, κ, cdl0)
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

# ╔═╡ a84ec6c7-6ee1-46e7-8a19-6a4c355fce6f
molvar_results = dlcap_calc(;
    nref = 1,
    molarities = [0.001, 0.01, 0.1, 1],
    κ = [10.0, 10, 10, 10]
)

# ╔═╡ f470297d-2649-424b-abe2-8418120a779a
dlcap_plot(molvar_results, "example_dlcap_molvar.png")

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
function bvcompare_calc(; nref = 0, zR = 1, n = 1, verbose = "", bgmol=2)
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
                boundary_dirichlet!(f, u, bnode; species = iϕ, region = Γ_we, value = ϕ_we)
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
    @info c_bulk' * z
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
        scheme
    )


    function sweep(; eneutral = true, tunnel = false, bikerman = true)
        celldata = deepcopy(pnpdata)
        celldata.eneutral = eneutral
        if !bikerman
            celldata.v .= 0.0
            celldata.κ .= 0.0
        end

        if tunnel
            pnpcell = PNPSystem(grid; bcondition = tunnelbc, reaction = tunnelreaction, celldata)
        else
            pnpcell = PNPSystem(grid; bcondition = halfcellbc, celldata)
        end
        result = ivsweep(pnpcell; voltages, store_solutions = true, verbose, Δu_opt = 1.0e10)
        currs = -zR * LiquidElectrolytes.currents(result, iR, electrode = :bulk) - zO * LiquidElectrolytes.currents(result, iO, electrode = :bulk)
        sol = LiquidElectrolytes.voltages_solutions(result)
        return result, currs, sol
    end

    pnpresult, pnpcurrs, pnpsol = sweep(eneutral = false, tunnel = false)
    nnpresult, nnpcurrs, nnpsol = sweep(eneutral = true, tunnel = false)
    return (; pnpresult, pnpcurrs, pnpsol, nnpresult, nnpcurrs, nnpsol, zR, n)
end


# ╔═╡ 21ab6d00-e979-4fce-bf6a-9cd4995d86f3
function bvcompare_plot(res, output = nothing; scaling = [1, 100, 1.0e5])
    (; pnpresult, pnpcurrs, pnpsol, nnpresult, nnpcurrs, nnpsol, zR, n) = res
    tunnel = false

    function yplot(fg, imax; xlabel = "", legend = nothing)
        ax = Axis(fg; ylabel = L"I/(mA/cm^2)", xlabel)
        ylims!(ax, [-imax, imax])
        lines!(
            ax,
            pnpresult.voltages, -pnpcurrs / (mA / cm^2), color = :red, label = "PNP"
        )
        # Plotter.scatter!(ax,pnpresult.voltages[1:12:end], -pnpcurrs[1:12:end]/(mA/cm^2),
        #          color = :red, markershape = :utriangle, markersize = 3)
        lines!(
            ax,
            nnpresult.voltages, -nnpcurrs / (mA / cm^2),
            color = :green, label = "NNP",
            linestyle = :dash,
            linewidth = 2
        )

        if tunnel
            lines!(p, tresult.voltages, -tcurrs / (mA / cm^2), color = :blue, label = "PNP-tunnel")
        end
        if !isnothing(legend)
            axislegend(;
                position = legend, labelsize = 10,
                backgroundcolor = RGBA(1.0, 1.0, 1.0, 0.5)
            )
        end
        return ax
    end
    fig = Figure(size = (600, 400), fontsize = 12)
    title = tunnel ? L"z_R=%$(zR), n=%$(n), β=%$(btxt)/cm" : L"z_R=%$(zR), n=%$(n)"
    yplot(fig[1, :], scaling[1], legend = :lt)
    yplot(fig[2, :], scaling[2])
    yplot(fig[3, :], scaling[3]; xlabel = L"U/V")
    Label(fig[0, :], title, tellwidth = false, fontsize = 16)

    if !isnothing(output)
        CairoMakie.save(joinpath(@__DIR__, output), fig)
    end
    return fig

end

# ╔═╡ 1d6728f1-c88f-48ab-85b4-afa9372bf1a3
bvc_resultm1 = bvcompare_calc(; zR = -1, n = 1)

# ╔═╡ 41a5b4ce-b793-4dd0-9a82-221abc12c5c5
bvcompare_plot(bvc_resultm1, "bvc_zr=-1.png", scaling = [1.5, 50, 1.0e5])

# ╔═╡ 45165906-d37a-4cec-8c98-2276f46054e3
bvc_result0 = bvcompare_calc(; zR = 0, n = 1)

# ╔═╡ 16539832-295d-49ef-a05c-be11d5d1b1a4
bvcompare_plot(bvc_result0, "bvc_zr=0.png", scaling = [1.5, 50, 1.0e5])

# ╔═╡ c58a3109-8eaf-4c55-86d4-f49158405b2f
bvc_resultp1 = bvcompare_calc(; zR = 1, n = 1)

# ╔═╡ 7080aff5-3f11-4c5d-8cfe-2f7995e95450
bvcompare_plot(bvc_resultp1, "bvc_zr=1.png", scaling = [0.5, 150, 1.0e5])

# ╔═╡ 955a53ac-5a2f-463f-bf10-ca50159861b7
bvc_resultx = bvcompare_calc(; zR = -1, n = 2, bgmol=0.1)

# ╔═╡ 920a3fc6-91d4-42f1-ad49-a446703b6c46
bvcompare_plot(bvc_resultx, "bvc_x.png", scaling = [0.5, 300, 5.0e5])

# ╔═╡ 8af12f1c-d35b-4cc9-8185-1bb5adbb69e8
html"""<hr>"""

# ╔═╡ 784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
html"""<style>.dont-panic{ display: none }</style>"""

# ╔═╡ Cell order:
# ╠═a70cef7d-2a2f-4155-bdf3-fec9df94c63f
# ╠═f674a3f5-87bf-4677-9081-4fd7b0cbb6e3
# ╟─b3f2b24f-2d3a-4095-94e7-1f858d7e1519
# ╠═2e6524ce-86f0-47cb-a344-c014638b0336
# ╠═b19bc908-1c88-4906-a9f7-3b240d9805b0
# ╠═e55f583a-a249-4cc1-b310-bfde0cab9d12
# ╠═a84ec6c7-6ee1-46e7-8a19-6a4c355fce6f
# ╠═f470297d-2649-424b-abe2-8418120a779a
# ╠═ad85052b-f5cd-451b-a978-7df924b1cd61
# ╠═17ba9584-d7a2-4745-94f7-5bd651e67fa5
# ╟─aef407a4-b3d4-403b-99b0-6fc29883ead9
# ╠═a1b85df4-954f-4c77-9d07-ab6e0f6d5745
# ╠═21ab6d00-e979-4fce-bf6a-9cd4995d86f3
# ╠═1d6728f1-c88f-48ab-85b4-afa9372bf1a3
# ╠═41a5b4ce-b793-4dd0-9a82-221abc12c5c5
# ╠═45165906-d37a-4cec-8c98-2276f46054e3
# ╠═16539832-295d-49ef-a05c-be11d5d1b1a4
# ╠═c58a3109-8eaf-4c55-86d4-f49158405b2f
# ╠═7080aff5-3f11-4c5d-8cfe-2f7995e95450
# ╠═955a53ac-5a2f-463f-bf10-ca50159861b7
# ╠═920a3fc6-91d4-42f1-ad49-a446703b6c46
# ╟─8af12f1c-d35b-4cc9-8185-1bb5adbb69e8
# ╟─784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
