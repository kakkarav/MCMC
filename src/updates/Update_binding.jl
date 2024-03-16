function metro_binding!(
    lat::Lattice,
    params::SimParams,
    disorder::Disorder,
    rng,
    location::NTuple{NDIMS,Int64},
    dir::Int64,
)
    """
    this is a special update that update P and current J simultaneously.
    by changing P, we create a loop in curl variable. 
    """
    #propose an update the current and the gauge field simultaneously 
    random_move = rand(rng, binding_range)
    p_change = mod1(random_move, 3) - 2
    loop_change = fld1(random_move, 3) - 2

    if (p_change != 0 || loop_change != 0)
        dE = energy_local_villain(lat, params, disorder, location, p_change, dir)
        dE +=
            energy_coupling_old(lat, params, disorder, location, p_change, loop_change, dir)
        dE += energy_loop_disorder(lat, params, disorder, location, loop_change, dir)
        #Metropolis transition probability
        if (dE <= 0.0 || rand(rng) <= exp.(-dE))
            lat.energy += dE
            lat.p_lattice[location..., dir] += p_change
            update_c!(lat, location, p_change, dir) #update curl variable
            update_binding!(lat, location, loop_change, dir) #update current loop
        end
    end
end