#update the fluctuating angle of Villain model
function metro_fluctuate!(
    lattice::Lattice,
    params::SimParams,
    disorder::Disorder,
    rng::Random.Xoshiro,
    dir::Int64,
)
    angle_change = 2.0 * (rand(rng) - 0.5) * params.delta_theta
    #There is no binding/ coupling energy related to the update of teh angle variable
    dE =
        Energy_local_fluctuate(lattice, params, disorder, dir, angle_change) -
        Energy_local_fluctuate(lattice, params, disorder, dir, 0.0)
    if (dE <= 0.0 || rand(rng) <= exp.(-dE))
        lattice.fluctuate_lattice[dir] =
            update_angle(lattice.fluctuate_lattice[dir], angle_change)
        lattice.energy += dE
    end
end

function Energy_local_fluctuate(
    lattice::Lattice,
    params::SimParams,
    disorder::Disorder,
    dir::Int64,
    delta_angle::Float64,
)
    E_tot = 0.0
    new_fluctuate = update_angle(lattice.fluctuate_lattice[dir], delta_angle) / params.L
    for site in eachindex(lattice.angle)
        location = Tuple(CartesianIndices(lattice.angle)[site])
        nn = hop(location, dir, 1, size(lattice.angle))
        bond =
            dir == 3 ? params.lambda1 :
            disorder.bond_disorder_villain[location[1], location[2], dir]
        E =
            bond *
            (
                lattice.angle[nn...] - lattice.angle[location...] -
                2 * pi * lattice.p_lattice[location..., dir] - new_fluctuate
            )^2.0 / 2.0
        E_tot += E
    end
    return E_tot
end