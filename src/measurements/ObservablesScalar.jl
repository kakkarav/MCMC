# Measure the number of closed looped or the partition function ratio
struct ObsZRatio <: Obs
    obs_data::ObsData{Float64}
    ObsZRatio(num_of_measure) = new(ObsDataScalar("ZRatio", num_of_measure))
end
function measure_obs!(obs_Zratio::ObsZRatio, sim)
    closed = sim.lat.head == sim.lat.tail ? 1.0 : 0.0
    take_measurement!(obs_Zratio.obs_data, closed)
end


###########################################################################################################
# Energy
struct ObsEnergy <: Obs
    obs_data::ObsData{Float64}
    ObsEnergy(num_of_measure) = new(ObsDataScalar("Energy", num_of_measure))
end

function measure_obs!(obs_E::ObsEnergy, sim)
    E = sim.lat.energy / sim.sim_params.L^(NDIMS - 1) / sim.sim_params.Lt
    take_measurement!(obs_E.obs_data, E)
    return nothing
end


################################################################################################################

#Compressibility: stiffnesss in the imaginary time direction
struct ObsCompressLoop <: Obs
    obs_data::ObsData{Float64}
    ObsCompressLoop(num_of_measure) =
        new(ObsDataScalar("Compressibility Loop", num_of_measure))
end

# TODO: find the right definition of compressibility in term of beta
function measure_obs!(obs::ObsCompressLoop, sim)
    par_numsq = (sim.lat.winding[end]) .^ 2.0
    compress = par_numsq / sim.sim_params.L^2 ./ sim.sim_params.Lt
    take_measurement!(obs.obs_data, compress)
    return nothing
end

###########################################################################################################
#stiffness loop
struct ObsStiffnessLoopAll <: Obs
    obs_data::ObsData{Array{Float64,1}}
    stiffFT::Array{Float64,1}
    current::Array{Float64,1}
    ObsStiffnessLoopAll(num_of_measure, Lt) = new(
        ObsDataArray("Stiffness Loop All", num_of_measure, ntuple(x -> Lt, 1)),
        zeros(ntuple(x -> Lt, 1)),
        zeros(ntuple(x -> Lt, 1)),
    )
end

function measure_obs!(obs::ObsStiffnessLoopAll, sim)
    fill!(obs.stiffFT, 0.0)
    for i = 1:sim.sim_params.Lt
        for dir = 1:2
            fill!(obs.current, 0.0)
            for t = 1:sim.sim_params.Lt
                obs.current[t] += cal_current_loop_list(sim, dir, t)
            end
            obs.stiffFT[i] +=
                abs2.(dot(obs.current, sim.fourier[i])) ./ sim.sim_params.L^(NDIMS - 1) ./
                sim.sim_params.Lt ./ 2
        end
    end
    take_measurement!(obs.obs_data, obs.stiffFT)
    return nothing
end

###########################################################################################################
#Measure current of loop model
struct ObsCurrentTimeLoop <: Obs
    obs_data::ObsData{Float64}
    ObsCurrentTimeLoop(num_of_measure) = new(ObsDataScalar("Current Time Loop", num_of_measure))
end

function measure_obs!(obs_current::ObsCurrentTimeLoop, sim)
    current = float(sim.lat.winding[3]) / sim.sim_params.L^((NDIMS - 1) / 2) / sim.sim_params.Lt^(1 / 2)
    take_measurement!(obs_current.obs_data, current)
    return nothing
end

###########################################################################################################
struct ObsCurrentSpaceLoop <: Obs
    obs_data::ObsData{Array{Float64,1}}
    green::Array{Float64,1}
    ObsCurrentSpaceLoop(num_of_measure) = new(
        ObsDataArray("Current Space Loop", num_of_measure, ntuple(x -> 2, 1)),
        zeros(ntuple(x -> 2, 1)),
    )
end

function measure_obs!(obs_green::ObsCurrentSpaceLoop, sim)
    for dir = 1:2
        obs_green.green[dir] =
            float(sim.lat.winding[dir]) / sim.sim_params.L^((NDIMS - 1) / 2) /
            sim.sim_params.Lt^(1 / 2)
    end
    take_measurement!(obs_green.obs_data, obs_green.green)
    return nothing
end

###########################################################################################################
#Measure current of loop model
struct ObsCurrentLoopFT <: Obs
    obs_data::ObsData{Complex{Float64}}
    current::Array{Float64,1}
    ObsCurrentLoopFT(num_of_measure, Lt) =
        new(ObsDataComplex("Current Loop FT", num_of_measure), zeros(Lt))
end

function measure_obs!(obs_current::ObsCurrentLoopFT, sim)
    fill!(obs_current.current, 0.0)
    for t = 1:sim.sim_params.Lt
        obs_current.current[t] += cal_current_loop_list(sim, 1, t)
    end
    current = dot(obs_current.current, sim.fourier[2])
    current = current / sim.sim_params.L^(NDIMS / 2)
    take_measurement!(obs_current.obs_data, current)
    return nothing
end

###########################################################################################################
#Zero-momentum loop stiffness
struct ObsStiffnessLoop <: Obs
    obs_data::ObsData{Float64}
    ObsStiffnessLoop(num_of_measure) = new(ObsDataScalar("Stiffness Loop", num_of_measure))
end

function measure_obs!(obs_stiffness::ObsStiffnessLoop, sim)
    stiffness = 0.0
    stiff_dir = 0.0
    dir = 1
    for site in eachindex(sim.lat.angle)
        location = Tuple(CartesianIndices(sim.lat.angle)[site])
        stiff_dir += sim.lat.current[location..., dir]
    end
    stiffness += stiff_dir^2 / sim.sim_params.L^(NDIMS - 1) ./ sim.sim_params.Lt
    take_measurement!(obs_stiffness.obs_data, stiffness)
    return nothing
end

###########################################################################################################
#fourier[2] transform version of ObsStiffnessLoop
struct ObsStiffnessLoopFT <: Obs
    obs_data::ObsData{Float64}
    current::Array{Float64,1}
    ObsStiffnessLoopFT(num_of_measure, Lt) =
        new(ObsDataScalar("Stiffness Loop FT", num_of_measure), zeros(Lt))
end

function measure_obs!(obs_stiff::ObsStiffnessLoopFT, sim)
    fill!(obs_stiff.current, 0.0)
    for t = 1:sim.sim_params.Lt
        obs_stiff.current[t] += cal_current_loop_list(sim, 1, t)
    end
    stiffness =
        abs2.(dot(obs_stiff.current, sim.fourier[2])) / sim.sim_params.L^(NDIMS - 1) ./
        sim.sim_params.Lt
    take_measurement!(obs_stiff.obs_data, stiffness)
    return nothing
end

################################################################################################################

#Current of villain_J
struct ObsCurrentVillainFT <: Obs
    obs_data::ObsData{Complex{Float64}}
    current_villain::Array{Float64,1}
    current_loop::Array{Float64,1}
    ObsCurrentVillainFT(num_of_measure, Lt) =
        new(ObsDataComplex("Current Villain FT", num_of_measure), zeros(Lt), zeros(Lt))
end

function measure_obs!(obs_current::ObsCurrentVillainFT, sim)
    fill!(obs_current.current_villain, 0.0)
    fill!(obs_current.current_loop, 0.0)
    for t = 1:sim.sim_params.Lt
        obs_current.current_villain[t] += cal_current_villain_list(sim, 1, t)
        obs_current.current_loop[t] += cal_current_disorder_loop_list(sim, 2, t)
    end
    current_vill = dot(obs_current.current_villain, sim.fourier[2])
    # current_loop = dot(obs_current.current_loop,sim.fourier[2])*2*im*sin(pi/sim.sim_params.L)*exp(-im*pi/sim.sim_params.L)/(2*pi)
    current_loop =
        dot(obs_current.current_loop, sim.fourier[2]) *
        2 *
        im *
        sin(pi / sim.sim_params.Lt) *
        sim.sim_params.eta / (2 * pi)
    current =
        (current_vill - current_loop) ./ sim.sim_params.L ./ sim.sim_params.Lt^(1.0 ./ 2.0)

    take_measurement!(obs_current.obs_data, current)
    return nothing
end

###########################################################################################################
struct ObsStiffnessVillain <: Obs
    obs_data::ObsData{Float64}
    ObsStiffnessVillain(num_of_measure) =
        new(ObsDataScalar("Stiffness Villain", num_of_measure))
end

function measure_obs!(obs_stiff::ObsStiffnessVillain, sim)
    stiffness = cal_stiffness_villain(sim)
    take_measurement!(obs_stiff.obs_data, stiffness)
    return nothing
end

function cal_stiffness_villain(sim)
    stiffness = 0.0
    for dir in 1:2
        sum_current = 0.0
        for site in eachindex(sim.lat.angle)
            location = Tuple(CartesianIndices(sim.lat.angle)[site])
            nn = hop(location, dir, 1, size(sim.lat.angle))
            bond = sim.disorder.bond_disorder_villain[location[1], location[2], dir]
            # fluctuate = location[dir]==1 ? sim.lat.fluctuate_lattice[dir] : 0.0
            fluctuate = sim.lat.fluctuate_lattice[dir] / sim.sim_params.L
            sum_current +=
                bond * (
                    sim.lat.angle[nn...] - sim.lat.angle[location...] -
                    2 * pi * sim.lat.p_lattice[location..., dir] - fluctuate
                )
        end
        stiffness +=
            (sim.disorder.sum_villain_bond[1] * sim.sim_params.Lt - sum_current^2) /
            sim.sim_params.L^(NDIMS - 1) ./ sim.sim_params.Lt
    end
    return stiffness / 2
end

###########################################################################################################
struct ObsCompressVillain <: Obs #compressibility
    obs_data::ObsData{Float64}
    ObsCompressVillain(num_of_measure) =
        new(ObsDataScalar("Compressibility Villain", num_of_measure))
end

function measure_obs!(obs_stiff::ObsCompressVillain, sim)
    compressibility = cal_compressibility_villain(sim)
    take_measurement!(obs_stiff.obs_data, compressibility)
    return nothing
end

function cal_compressibility_villain(sim)
    stiffness = 0.0
    sum_current = 0.0
    dir = 3
    for site in eachindex(sim.lat.angle)
        location = Tuple(CartesianIndices(sim.lat.angle)[site])
        # nn = hop(location,dir,1,size(sim.lat.angle))
        bond = sim.sim_params.lambda1
        # fluctuate = location[dir]==1 ? sim.lat.fluctuate_lattice[dir] : 0.0
        fluctuate = sim.lat.fluctuate_lattice[dir] / sim.sim_params.L
        # sum_current += bond*(sim.lat.angle[nn...]-sim.lat.angle[location...]-2*pi*sim.lat.p_lattice[location...,dir]-fluctuate)
        sum_current += bond * (-2 * pi * sim.lat.p_lattice[location..., dir] - fluctuate)
    end
    stiffness +=
        (
            sim.sim_params.lambda1 * (sim.sim_params.L)^2 * sim.sim_params.Lt -
            sum_current^2
        ) ./ sim.sim_params.L^(NDIMS - 1) ./ sim.sim_params.Lt
    return stiffness
end

###########################################################################################################
#calculate the current in the time direction
struct ObsCurrentTimeVillain <: Obs
    obs_data::ObsData{Float64}
    ObsCurrentTimeVillain(num_of_measure) =
        new(ObsDataScalar("Current Time Villain", num_of_measure))
end

function measure_obs!(obs_stiff::ObsCurrentTimeVillain, sim)
    current = cal_current_time_villain(sim)
    take_measurement!(obs_stiff.obs_data, current)
    return nothing
end

function cal_current_time_villain(sim)
    current = 0.0
    sum_current = 0.0
    dir = 3
    for site in eachindex(sim.lat.angle)
        location = Tuple(CartesianIndices(sim.lat.angle)[site])
        nn = hop(location, dir, 1, size(sim.lat.angle))
        bond = sim.sim_params.lambda1
        sum_current +=
            bond * (
                sim.lat.angle[nn...] - sim.lat.angle[location...] -
                2 * pi * sim.lat.p_lattice[location..., dir]
            )
    end
    current += sum_current
    return current ./ sim.sim_params.L ./ sim.sim_params.Lt^(1 / 2)
end


###########################################################################################################

#Fourier transform version of ObsStiffnessVillain
struct ObsStiffnessVillainFT <: Obs
    obs_data::ObsData{Float64}
    current_villain::Array{Float64,1}
    current_loop::Array{Float64,1}
    ObsStiffnessVillainFT(num_of_measure, Lt) =
        new(ObsDataScalar("Stiffness Villain FT", num_of_measure), zeros(Lt), zeros(Lt))
end

function measure_obs!(obs_stiff::ObsStiffnessVillainFT, sim)
    fill!(obs_stiff.current_villain, 0.0)
    fill!(obs_stiff.current_loop, 0.0)
    for t = 1:sim.sim_params.Lt
        obs_stiff.current_villain[t] += cal_current_villain_list(sim, 1, t)
        obs_stiff.current_loop[t] += cal_current_disorder_loop_list(sim, 2, t)
    end
    current_vill = dot(obs_stiff.current_villain, sim.fourier[2])
    # current_loop = dot(obs_stiff.current_loop,sim.fourier)*2*im*sin(pi/sim.sim_params.L)*exp(-im*pi/sim.sim_params.L)/(2*pi)
    current_loop =
        dot(obs_stiff.current_loop, sim.fourier[2]) *
        2 *
        im *
        sin(pi / sim.sim_params.Lt) *
        sim.sim_params.eta / (2 * pi)

    correlation = abs2.(current_vill - current_loop)

    stiffness =
        (
            (
                sim.disorder.sum_villain_bond[1] * sim.sim_params.Lt +
                4 *
                sin(pi / sim.sim_params.Lt)^2.0 *
                sim.disorder.sum_loop_bond_inverse[2] *
                sim.sim_params.Lt / (2 * pi)^2
            ) - correlation
        ) ./ sim.sim_params.L^(NDIMS - 1) ./ sim.sim_params.Lt
    take_measurement!(obs_stiff.obs_data, stiffness)
    return nothing
end

###########################################################################################################
struct ObsMagnitization <: Obs
    obs_data::ObsData{Float64}
    ObsMagnitization(num_of_measure) = new(ObsDataScalar("Magnitization", num_of_measure))
end

function measure_obs!(obs_mag::ObsMagnitization, sim)
    tot_mag = 0.0 + im * 0.0
    for site in eachindex(sim.lat.angle)
        location = Tuple(CartesianIndices(sim.lat.angle)[site])
        tot_mag += exp.(im * sim.lat.angle[location...])
    end
    magnitization = abs.(tot_mag) ./ sim.sim_params.L^(NDIMS - 1) ./ sim.sim_params.Lt
    take_measurement!(obs_mag.obs_data, magnitization)
    return nothing
end
