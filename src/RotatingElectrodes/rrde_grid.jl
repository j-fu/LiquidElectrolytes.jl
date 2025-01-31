"""
Grids for RRDE discretization
"""


const b_inert = 1  # inert part of disk at bottom and gap
const b_in = 2     # inlet at top
const b_out = 3    # outflow  through cylinder mantle
const b_disk = 4   # reactive part of disk
const b_ring = 5   # ring electrode
const b_symm = 6   # symmetry (r=0)
const b_max = b_symm

function rrde_grid(
        geom;
        nref = 0,
        hz = (0.05, 0.05),
        esmall = 0.015,
        elarge = 0.15,
        exlarge = 0.4
    )

    allmaskmin = [0.0, 0.0, 0.0]
    allmaskmax = [geom.r_4, geom.h_cyl, 2.0 * 3.15]

    symmmaskmin = [0.0, 0.0, 0.0]
    symmmaskmax = [0.0, geom.h_cyl, 0.0]

    inmaskmin = [0.0, geom.h_cyl, 0.0]
    inmaskmax = [geom.r_4, geom.h_cyl, 2.0 * 3.15]

    diskmaskmin = [0.0, 0.0, 0.0]
    diskmaskmax = [geom.r_1, 0.0, 2.0 * 3.15]

    ringmaskmin = [geom.r_2, 0.0, 0.0]
    ringmaskmax = [geom.r_3, 0.0, 2.0 * 3.15]

    outmaskmin = [geom.r_4, 0, 0.0]
    outmaskmax = [geom.r_4, geom.h_cyl, 2 * 3.15]

    D = 1.0 / 2.0^(nref)
    d_ring = (geom.r_3 - geom.r_2)
    d_gap = (geom.r_2 - geom.r_1)
    d_shr = (geom.r_4 - geom.r_3)

    linspace(a, b, n) = collect(a:((b - a) / n):b)

    rcoord = geomspace(0.0, geom.r_1, D * exlarge * geom.r_1, D * esmall * d_gap)
    startPawel = length(rcoord)
    rcoord = glue(rcoord, geomspace(geom.r_1, 0.5 * (geom.r_1 + geom.r_2), D * esmall * d_gap, D * elarge * d_gap))
    rcoord = glue(rcoord, geomspace(0.5 * (geom.r_1 + geom.r_2), geom.r_2, D * elarge * d_gap, D * esmall * d_gap))
    rcoord = glue(rcoord, geomspace(geom.r_2, 0.5 * (geom.r_2 + geom.r_3), D * esmall * d_ring, D * elarge * d_ring))
    rcoord = glue(rcoord, geomspace(0.5 * (geom.r_2 + geom.r_3), geom.r_3, D * elarge * d_ring, D * esmall * d_ring))
    endPawel = length(rcoord)
    rcoord = glue(rcoord, geomspace(geom.r_3, geom.r_4 - geom.w_gap, D * esmall * d_ring, D * exlarge * d_shr))
    if geom.w_gap > 0
        rcoord = glue(rcoord, linspace(geom.r_4 - geom.w_gap, geom.r_4, 2))
    end
    rSize = length(rcoord)
    zcoord = geomspace(0.0, geom.h_cyl, D * hz[1] * geom.h_cyl, D * hz[2] * geom.h_cyl)

    println("Points inbetween gap and ring: ", endPawel - startPawel)
    println("Points in r-direction: $rSize, i.e. Boxes in r-Direction: ", rSize - 1)
    println("Points in z-direction: ", length(zcoord))


    grid = VoronoiFVM.Grid(rcoord, zcoord)
    bfacemask!(grid, allmaskmin, allmaskmax, b_inert, tol = 1.0e-10, allow_new = false)
    bfacemask!(grid, inmaskmin, inmaskmax, b_in, tol = 1.0e-10, allow_new = false)
    bfacemask!(grid, outmaskmin, outmaskmax, b_out, tol = 1.0e-10, allow_new = false)
    bfacemask!(grid, diskmaskmin, diskmaskmax, b_disk, tol = 1.0e-10, allow_new = false)
    bfacemask!(grid, ringmaskmin, ringmaskmax, b_ring, tol = 1.0e-10, allow_new = false)
    bfacemask!(grid, symmmaskmin, symmmaskmax, b_symm, tol = 1.0e-10, allow_new = false)

    circular_symmetric!(grid)
    return grid
end


function rescale_z!(grid, delta)
    coord = coordinates(grid)
    maxz = coord[2, end]
    new_zscale = delta * 2 / maxz
    for i in 1:size(coord, 2)
        coord[2, i] = coord[2, i] * new_zscale
    end
    return grid[YCoordinates] .*= new_zscale
end


function rrde_grid_zgeom(geom, n_ref)
    allmaskmin = [0.0, 0.0, 0.0]
    allmaskmax = [geom.r_4, geom.h_cyl, 2.0 * 3.15]

    symmmaskmin = [0.0, 0.0, 0.0]
    symmmaskmax = [0.0, geom.h_cyl, 0.0]

    inmaskmin = [0.0, geom.h_cyl, 0.0]
    inmaskmax = [geom.r_4 - geom.w_gap, geom.h_cyl, 2.0 * 3.15]

    diskmaskmin = [0.0, 0.0, 0.0]
    diskmaskmax = [geom.r_1, 0.0, 2.0 * 3.15]

    ringmaskmin = [geom.r_2, 0.0, 0.0]
    ringmaskmax = [geom.r_3, 0.0, 2.0 * 3.15]

    outmaskmin = [geom.r_4, 0, 0.0]
    outmaskmax = [geom.r_4, geom.h_cyl, 2 * 3.15]

    D = 1.0 / 2.0^(n_ref)
    d_ring = (geom.r_3 - geom.r_2)
    d_gap = (geom.r_2 - geom.r_1)
    d_shr = (geom.r_4 - geom.r_3)
    esmall = 0.025
    elarge = 0.2
    exlarge = 0.4
    linspace(a, b, n) = collect(a:((b - a) / n):b)

    rcoord = geomspace(0.0, geom.r_1, D * exlarge * geom.r_1, D * esmall * d_gap)
    rcoord = glue(rcoord, geomspace(geom.r_1, 0.5 * (geom.r_1 + geom.r_2), D * esmall * d_gap, D * elarge * d_gap))
    rcoord = glue(rcoord, geomspace(0.5 * (geom.r_1 + geom.r_2), geom.r_2, D * elarge * d_gap, D * esmall * d_ring))
    rcoord = glue(rcoord, geomspace(geom.r_2, 0.5 * (geom.r_2 + geom.r_3), D * esmall * d_ring, D * elarge * d_ring))
    rcoord = glue(rcoord, geomspace(0.5 * (geom.r_2 + geom.r_3), geom.r_3, D * elarge * d_ring, D * esmall * d_ring))
    rcoord = glue(rcoord, geomspace(geom.r_3, geom.r_4 - geom.w_gap, D * esmall * d_ring, D * exlarge * d_shr))

    if geom.w_gap > 0
        rcoord = glue(rcoord, linspace(geom.r_4 - geom.w_gap, geom.r_4, 2))
    end

    zcoord = geomspace(0.0, geom.h_cyl, D * 0.01 * geom.h_cyl, D * 0.25 * geom.h_cyl)

    grid = VoronoiFVM.Grid(rcoord, zcoord)
    bfacemask!(grid, allmaskmin, allmaskmax, b_inert, tol = 1.0e-10)
    bfacemask!(grid, inmaskmin, inmaskmax, b_in, tol = 1.0e-10)
    bfacemask!(grid, outmaskmin, outmaskmax, b_out, tol = 1.0e-10)
    bfacemask!(grid, diskmaskmin, diskmaskmax, b_disk, tol = 1.0e-10)
    bfacemask!(grid, ringmaskmin, ringmaskmax, b_ring, tol = 1.0e-10)
    bfacemask!(grid, symmmaskmin, symmmaskmax, b_symm, tol = 1.0e-10)

    circular_symmetric!(grid)
    return grid
end
