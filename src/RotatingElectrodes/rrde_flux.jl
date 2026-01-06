mutable struct RRDEData
    nspecies::Int64
    D::Vector{Float64}
    bm::Matrix{Float64}
    bp::Matrix{Float64}
    bfvelo::Matrix{Float64}
    RRDEData(diffcoeffs::AbstractVector) = RRDEData(new(), diffcoeffs)
end

function RRDEData(data, diffcoeffs::AbstractVector)
    data.nspecies = length(diffcoeffs)
    data.D = convert(Vector{Float64}, copy(diffcoeffs))
    return data
end

function update!(data::RRDEData, grid, func)
    evelo = edgevelocities(grid, func)
    data.bfvelo = bfacevelocities(grid, func)
    data.bm = zeros(data.nspecies, length(evelo))
    data.bp = zeros(data.nspecies, length(evelo))
    for idx in 1:length(evelo)
        for ispec in 1:data.nspecies
            xbp, xbm = VoronoiFVM.fbernoulli_pm(evelo[idx] / data.D[ispec])
            data.bp[ispec, idx] = data.D[ispec] * xbp
            data.bm[ispec, idx] = data.D[ispec] * xbm
        end
    end
    return
end

function rrde_flux(f, u, edge, data)
    idx = edge.index
    for ispec in 1:data.nspecies
        f[ispec] = data.bp[ispec, idx] * u[ispec, 1] - data.bm[ispec, idx] * u[ispec, 2]
    end
    return
end

function rrde_outflow(f, u, node, data)
    return if node.region == b_out
        v = data.bfvelo[node.ibnode, node.ibface]
        for ispec in 1:data.nspecies
            f[ispec] = v * u[ispec]
        end
    end
end
