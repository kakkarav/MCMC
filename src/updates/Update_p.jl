function metro_p!(
    lat::Lattice,
    params::SimParams,
    disorder::Disorder,
    rng,
    location::NTuple{NDIMS,Int64},
    dir::Int64,
)
    #change p by +-1
    p_change = 2 * rand(rng, 0:1) - 1
    dE_1 = energy_local_villain(lat, params, disorder, location, p_change, dir)
    dE_2 = energy_coupling_old(lat, params, disorder, location, p_change, 0, dir)
    dE = dE_1 + dE_2
    if (dE <= 0.0 || rand(rng) <= exp.(-dE))
        lat.energy += dE
        lat.p_lattice[location..., dir] += p_change
        update_c!(lat, location, p_change, dir)
    end
    return nothing
end

function energy_local_villain(
    lat::Lattice,
    params::SimParams,
    disorder::Disorder,
    location::NTuple{NDIMS,Int64},
    p_change::Int64,
    dir::Int64,
)
    """
    Compute Energy change from the Villain part (no binding)
    """
    nn = hop(location, dir, 1, size(lat.angle))
    angle = lat.angle
    # fluctuate = location[dir]==1 ? lat.fluctuate_lattice[dir] : 0.0
    fluctuate = lat.fluctuate_lattice[dir] / params.L
    p_lattice = lat.p_lattice
    lambda1 =
        dir == 3 ? params.lambda1 :
        disorder.bond_disorder_villain[location[1], location[2], dir]
    #only one bond is altered by changing P, so there is only one bond changed for Villain part
    E_new =
        (
            angle[nn...] - angle[location...] -
            2 * pi * (p_lattice[location..., dir] + p_change) - fluctuate
        )^2.0
    E_old =
        (
            angle[nn...] - angle[location...] - 2 * pi * p_lattice[location..., dir] -
            fluctuate
        )^2.0
    return (E_new - E_old) * lambda1 / 2.0
end