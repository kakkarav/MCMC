#this will calculate the energy of a plaquette in the curl and current variable
#this shows up in update_binding because flipping p create a loop
function energy_coupling_old(
    lat::Lattice,
    sim_params::SimParams,
    disorder::Disorder,
    location::NTuple{NDIMS,Int64},
    p_change::Int64,
    loop_change::Int64,
    d::Int64,
)
    #The update I did below here is a bit obscure, but can be explain by drawing a picture. Here I assume that the relative vector we between the p_lattice (villain model) and the dual lattice(loop model) is
    # r_loop - r_villain  = (+x,+y,+z), so I pick the dual lattice to be offset from the villain lattice in the positive direction.
    #Below I did two jumps because of the above definition of relative direction of the two lattices
    #location is in the Villain coordinate

    d_1 = d                # z
    d_2 = mod1(d + 1, NDIMS)   # x
    d_3 = mod1(d + 2, NDIMS)   # y

    loop_location = location
    eta = sim_params.eta
    E = 0.0

    if (d == 3)
        #no chemical potential in this case
        # #change location into the loop coordinate
        # loop_location1 = hop(location, d_2, 2, size(lat.angle))
        # loop_location = hop(loop_location1, d_3, 2, size(lat.angle))

        #\delta curl has different sign for each site
        # jump to y -- > y-1
        jump_1 = hop(loop_location, d_3, 2, size(lat.angle))
        bond_1 = disorder.bond_disorder_loop[jump_1[1], jump_1[2], d_3]

        # compute energy diff in the y direction -- increase charge
        E +=
            (
                (loop_change - eta * p_change)^2 +
                2 *
                (loop_change - eta * p_change) *
                (lat.current[jump_1..., d_3] - eta * lat.curl_lattice[jump_1..., d_3])
            ) / (2.0 * bond_1)

        # jump to x -- > x-1
        jump_2 = hop(loop_location, d_2, 2, size(lat.angle))
        bond_2 = disorder.bond_disorder_loop[jump_2[1], jump_2[2], d_2]
        # compute energy diff in the x direction -- decrease charge
        E +=
            (
                (loop_change - eta * p_change)^2 -
                2 *
                (loop_change - eta * p_change) *
                (lat.current[jump_2..., d_2] - eta * lat.curl_lattice[jump_2..., d_2])
            ) / (2.0 * bond_2)

        # jump to x -- > x-1 and y --> y-1
        jump_3 = hop(jump_2, d_3, 2, size(lat.angle))
        bond_3 = disorder.bond_disorder_loop[jump_3[1], jump_3[2], d_2]
        bond_4 = disorder.bond_disorder_loop[jump_3[1], jump_3[2], d_3]
        # compute energy diff in the x direction -- increase charge
        # compute energy diff in the y direction -- decrease charge

        E +=
            (
                (loop_change - eta * p_change)^2 +
                2 *
                (loop_change - eta * p_change) *
                (lat.current[jump_3..., d_2] - eta * lat.curl_lattice[jump_3..., d_2])
            ) / (2.0 * bond_3)
        E +=
            (
                (loop_change - eta * p_change)^2 -
                2 *
                (loop_change - eta * p_change) *
                (lat.current[jump_3..., d_3] - eta * lat.curl_lattice[jump_3..., d_3])
            ) / (2.0 * bond_4)

    elseif (d == 1)
        #d=x,d_2=y,d_3=z

        #\delta curl has different sign for each site
        # jump to z -- > z-1
        jump_1 = hop(loop_location, 3, 2, size(lat.angle))
        bond_1 = sim_params.lambda2
        # compute energy diff in the y direction -- increase charge
        E +=
            (
                (eta * p_change - loop_change)^2 +
                2 *
                (eta * p_change - loop_change) *
                (eta * lat.curl_lattice[jump_1..., 3] - lat.current[jump_1..., 3])
            ) / (2.0 * bond_1)

        # jump to y -- > y-1
        jump_2 = hop(loop_location, 2, 2, size(lat.angle))
        bond_2 = disorder.bond_disorder_loop[jump_2[1], jump_2[2], 2]
        # compute energy diff in the y direction -- decrease charge
        E +=
            (
                (eta * p_change - loop_change)^2 -
                2 *
                (eta * p_change - loop_change) *
                (eta * lat.curl_lattice[jump_2..., 2] - lat.current[jump_2..., 2])
            ) / (2.0 * bond_2)

        # jump to z -- > z-1 and y --> y-1
        jump_3 = hop(jump_2, 3, 2, size(lat.angle))
        bond_3 = disorder.bond_disorder_loop[jump_3[1], jump_3[2], 2]
        bond_4 = sim_params.lambda2
        # compute energy diff in the y direction -- increase charge
        # compute energy diff in the z direction -- decrease charge

        E +=
            (
                (eta * p_change - loop_change)^2 +
                2 *
                (eta * p_change - loop_change) *
                (eta * lat.curl_lattice[jump_3..., 2] - lat.current[jump_3..., 2])
            ) / (2.0 * bond_3)
        E +=
            (
                (eta * p_change - loop_change)^2 -
                2 *
                (eta * p_change - loop_change) *
                (eta * lat.curl_lattice[jump_3..., 3] - lat.current[jump_3..., 3])
            ) / (2.0 * bond_4)


    elseif (d == 2)
        #d=y,d_2=z,d_3=x

        #\delta curl has different sign for each site
        # jump to z -- > z-1
        jump_1 = hop(loop_location, 3, 2, size(lat.angle))
        bond_1 = sim_params.lambda2
        # compute energy diff in the y direction -- decrease charge
        E +=
            (
                (eta * p_change - loop_change)^2 -
                2 *
                (eta * p_change - loop_change) *
                (eta * lat.curl_lattice[jump_1..., 3] - lat.current[jump_1..., 3])
            ) / (2.0 * bond_1)

        # jump to x -- > x-1
        jump_2 = hop(loop_location, 1, 2, size(lat.angle))
        bond_2 = disorder.bond_disorder_loop[jump_2[1], jump_2[2], 1]
        # compute energy diff in the x direction -- increase charge
        E +=
            (
                (eta * p_change - loop_change)^2 +
                2 *
                (eta * p_change - loop_change) *
                (eta * lat.curl_lattice[jump_2..., 1] - lat.current[jump_2..., 1])
            ) / (2.0 * bond_2)

        # jump to z -- > z-1 and x --> x-1
        jump_3 = hop(jump_2, 3, 2, size(lat.angle))
        bond_3 = disorder.bond_disorder_loop[jump_3[1], jump_3[2], 1]
        bond_4 = sim_params.lambda2
        # compute energy diff in the x direction -- decrease charge
        # compute energy diff in the z direction -- increase charge

        E +=
            (
                (eta * p_change - loop_change)^2 -
                2 *
                (eta * p_change - loop_change) *
                (eta * lat.curl_lattice[jump_3..., 1] - lat.current[jump_3..., 1])
            ) / (2.0 * bond_3)
        E +=
            (
                (eta * p_change - loop_change)^2 +
                2 *
                (eta * p_change - loop_change) *
                (eta * lat.curl_lattice[jump_3..., 3] - lat.current[jump_3..., 3])
            ) / (2.0 * bond_4)

    end

    return float(E)
end


#the function calculate the energy from the disordered chemical potential
function energy_loop_disorder(
    lat::Lattice,
    sim_params::SimParams,
    disorder::Disorder,
    location::NTuple{NDIMS,Int64},
    loop_change::Int64,
    d::Int64,
)
    E = 0.0
    charge = loop_change
    if (d == 1) #if a loop is not in the XY plane
        #if the loop is in the XY plane, there is no energy change from chemical potential
        #Else, a loop will have two counterpropagating currents in the time direction. Each of them is at different locations and experience different disorders.

        #current +1 for current at locations
        #current -1 for current at jump_1, y -> y-1
        jump_1 = hop(location, 3, 2, size(lat.angle))
        jump_2 = hop(jump_1, 2, 2, size(lat.angle))

        dE1 =
            charge * disorder.chemical_potential[jump_1[1], jump_1[2]] ./ sim_params.lambda2
        dE2 =
            (-charge) * disorder.chemical_potential[jump_2[1], jump_2[2]] ./
            sim_params.lambda2
        E += dE1 + dE2

    elseif (d == 2) # a loop is in the ZX plane
        jump_1 = hop(location, 3, 2, size(lat.angle))
        jump_2 = hop(jump_1, 1, 2, size(lat.angle))

        dE1 =
            (-charge) * disorder.chemical_potential[jump_1[1], jump_1[2]] ./
            sim_params.lambda2
        dE2 =
            charge * disorder.chemical_potential[jump_2[1], jump_2[2]] ./ sim_params.lambda2
        E += dE1 + dE2
    end
    return E
end



#update curl lattice
#change in P creates a loop in variable curl(P)
function update_c!(lat::Lattice, location::NTuple{NDIMS,Int64}, p_change::Int64, d::Int64)
    d_1 = d                # x
    d_2 = mod1(d + 1, NDIMS)   # y
    d_3 = mod1(d + 2, NDIMS)   # z

    # jump to z -- > z-1
    jump_1 = hop(location, d_3, 2, size(lat.angle))
    # compute energy diff in the z direction -- increase charge
    lat.curl_lattice[jump_1..., d_3] += p_change

    # jump to y -- > y-1
    jump_2 = hop(location, d_2, 2, size(lat.angle))
    # compute energy diff in the y direction -- decrease charge
    lat.curl_lattice[jump_2..., d_2] -= p_change

    # jump to y -- > y-1 and z --> z-1
    jump_3 = hop(jump_2, d_3, 2, size(lat.angle))
    # compute energy diff in the y direction -- increase charge
    # compute energy diff in the z direction -- decrease charge

    lat.curl_lattice[jump_3..., d_2] += p_change
    lat.curl_lattice[jump_3..., d_3] -= p_change

    # print(lat.curl_lattice[jump_1...,d_3],lat.curl_lattice[jump_2...,d_2],lat.curl_lattice[jump_3...,d_2],lat.curl_lattice[jump_3...,d_3])
    return nothing
end

#this update create a loop in J variable
function update_binding!(
    lat::Lattice,
    location::NTuple{NDIMS,Int64},
    loop_change::Int64,
    d::Int64,
)
    d_1 = d                # x
    d_2 = mod1(d + 1, NDIMS)   # y
    d_3 = mod1(d + 2, NDIMS)   # z

    # jump to z -- > z-1
    jump_1 = hop(location, d_3, 2, size(lat.angle))
    # compute energy diff in the z direction -- increase charge
    lat.current[jump_1..., d_3] += loop_change

    # jump to y -- > y-1
    jump_2 = hop(location, d_2, 2, size(lat.angle))
    # compute energy diff in the y direction -- decrease charge
    lat.current[jump_2..., d_2] -= loop_change

    # jump to y -- > y-1 and z --> z-1
    jump_3 = hop(jump_2, d_3, 2, size(lat.angle))
    # compute energy diff in the y direction -- increase charge
    # compute energy diff in the z direction -- decrease charge

    lat.current[jump_3..., d_2] += loop_change
    lat.current[jump_3..., d_3] -= loop_change

    return nothing
end