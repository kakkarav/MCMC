# Hall conductivity/ Cross specie stiffness
struct ObsHallConductivity <: Obs
  obs_data::ObsData{Complex{Float64}}
  villain_J::Array{Float64,1}
  loop_J::Array{Float64,1}
  coupled_J::Array{Float64,1}
  ObsHallConductivity(num_of_measure, Lt) = new(
    ObsDataComplex("Complex Hall Conductivity", num_of_measure),
    zeros(Lt),
    zeros(Lt),
    zeros(Lt),
  )
end


function measure_obs!(obs_hall_cond::ObsHallConductivity, sim)
  fill!(obs_hall_cond.villain_J, 0.0)
  fill!(obs_hall_cond.loop_J, 0.0)
  fill!(obs_hall_cond.coupled_J, 0.0)
  for t = 1:sim.sim_params.Lt
    obs_hall_cond.villain_J[t] += cal_current_villain_list(sim, 1, t)
    obs_hall_cond.loop_J[t] += cal_current_loop_list(sim, 2, t)
    obs_hall_cond.coupled_J[t] += cal_current_disorder_loop_list(sim, 2, t)
  end

  current_vill = dot(obs_hall_cond.villain_J, sim.fourier[2])
  current_coupled =
    dot(obs_hall_cond.coupled_J, sim.fourier[2]) *
    2 *
    im *
    sin(pi / sim.sim_params.Lt) *
    exp(-im * pi / sim.sim_params.Lt) *
    sim.sim_params.eta / (2 * pi)
  current_loop = conj(dot(obs_hall_cond.loop_J, sim.fourier[2]))
  correlation =
    (current_vill - current_coupled) * current_loop / sim.sim_params.L^(NDIMS - 1) /
    sim.sim_params.Lt
  correlation = correlation * exp(im * pi / sim.sim_params.Lt)

  take_measurement!(obs_hall_cond.obs_data, correlation)
  return nothing
end

#calculate a sum of the current of one time slice (t = time) for Villain model
function cal_current_villain_list(sim, dir::Int64, time::Int64)
  current = 0.0
  for i in 1:size(sim.lat.angle, 1)
    for j in 1:size(sim.lat.angle, 2)
      location = (i, j, time)
      nn = hop(location, dir, 1, size(sim.lat.angle))
      fluctuate = sim.lat.fluctuate_lattice[dir] / sim.sim_params.Lt
      current +=
        sim.disorder.bond_disorder_villain[i, j, dir] * (
          sim.lat.angle[nn...] - sim.lat.angle[location...] -
          2 .* pi * sim.lat.p_lattice[location..., dir] - fluctuate
        )
    end
  end
  return current
end

function cal_angle_villain_list(sim, dir::Int64, time::Int64)
  angle = 0.0
  for i in 1:size(sim.lat.angle, 1)
    for j in 1:size(sim.lat.angle, 2)
      location = (i, j, time)
      angle +=
        sim.disorder.bond_disorder_villain[i, j, dir] * sim.lat.angle[location...]
    end
  end
  return angle
end

function cal_current_curl_list(sim, dir::Int64, time::Int64)
  current = 0.0
  for i in 1:size(sim.lat.angle, 1)
    for j in 1:size(sim.lat.angle, 2)
      current += sim.lat.curl_lattice[i, j, time, dir]
    end
  end
  return current
end

#calculate a sum of the current of one time slice (t = time) for loop model
function cal_current_loop_list(sim, dir::Int64, time::Int64)
  current = 0.0
  for i in 1:size(sim.lat.angle, 1)
    for j in 1:size(sim.lat.angle, 2)
      current += sim.lat.current[i, j, time, dir]
    end
  end
  return current
end

function cal_current_disorder_loop_list(sim, dir::Int64, time::Int64)
  current = 0.0
  for i in 1:size(sim.lat.angle, 1)
    for j in 1:size(sim.lat.angle, 2)
      current +=
        (
          sim.lat.current[i, j, time, dir] -
          sim.sim_params.eta * sim.lat.curl_lattice[i, j, time, dir]
        ) / sim.disorder.bond_disorder_loop[i, j, dir]
    end
  end
  return current
end

###########################################################################################################

#single particle Green's function for loop model
struct ObsGreenLoopTime <: Obs
  obs_data::ObsData{Array{Float64,1}}
  green::Array{Float64,1}
  ObsGreenLoopTime(num_of_measure, Lt) = new(
    ObsDataArray("Green Loop Time", num_of_measure, ntuple(x -> Lt, 1)),
    zeros(ntuple(x -> Lt, 1)),
  )
end

# TODO: we need to modify so that we can find C(0,t) and C(r,0) separately
function measure_obs!(obs_green::ObsGreenLoopTime, sim)
  fill!(obs_green.green, 0.0)
  #we shift the index by 1 because julia start counting from 1 not 0
  delta_worm = mod(sim.lat.head[3] - sim.lat.tail[3], sim.sim_params.Lt) + 1
  if (sim.lat.head[1] == sim.lat.tail[1]) && (sim.lat.head[2] == sim.lat.tail[2])
    obs_green.green[delta_worm] += 1.0
  end
  # println(obs_green.green)
  take_measurement!(obs_green.obs_data, obs_green.green)
  return nothing
end

###########################################################################################################

struct ObsGreenLoopSpace <: Obs
  obs_data::ObsData{Array{Float64,1}}
  green::Array{Float64,1}
  ObsGreenLoopSpace(num_of_measure, L) = new(
    ObsDataArray("Green Loop Space", num_of_measure, ntuple(x -> L, 1)),
    zeros(ntuple(x -> L, 1)),
  )
end

# TODO: we need to modify so that we can find C(0,t) and C(r,0) separately
function measure_obs!(obs_green::ObsGreenLoopSpace, sim)
  fill!(obs_green.green, 0.0)
  #we shift the index by 1 because julia start counting from 1 not 0
  if (sim.lat.head[1] == sim.lat.tail[1])
    if (sim.lat.head[2] == sim.lat.tail[2])
      delta_worm = mod(sim.lat.head[1] - sim.lat.tail[1], sim.sim_params.Lt) + 1
      obs_green.green[delta_worm] += 1.0
    elseif (sim.lat.head[2] == sim.lat.tail[2])
      delta_worm = mod(sim.lat.head[1] - sim.lat.tail[1], sim.sim_params.Lt) + 1
      obs_green.green[delta_worm] += 1.0
    end
  end
  take_measurement!(obs_green.obs_data, obs_green.green)
  return nothing
end

###########################################################################################################

struct ObsCurrentSpaceVillain <: Obs
  obs_data::ObsData{Array{Float64,1}}
  green::Array{Float64,1}
  ObsCurrentSpaceVillain(num_of_measure) = new(
    ObsDataArray("Current Space Villain", num_of_measure, ntuple(x -> 2, 1)),
    zeros(ntuple(x -> 2, 1)),
  )
end

function measure_obs!(obs_green::ObsCurrentSpaceVillain, sim)
  for dir in 1:2
    obs_green.green[dir] = cal_current_villain(sim, dir)
  end
  take_measurement!(obs_green.obs_data, obs_green.green)
  return nothing
end

function cal_current_villain(sim, dir)
  current = 0.0
  sum_current = 0.0
  for site in eachindex(sim.lat.angle)
    location = Tuple(CartesianIndices(sim.lat.angle)[site])
    nn = hop(location, dir, 1, size(sim.lat.angle))
    bond = sim.disorder.bond_disorder_villain[location[1], location[2], dir]
    sum_current +=
      bond * (
        sim.lat.angle[nn...] - sim.lat.angle[location...] -
        2 * pi * sim.lat.p_lattice[location..., dir]
      )
  end
  current += sum_current
  return current / sim.sim_params.L^(2.0 ./ 2) ./ sim.sim_params.Lt^(1.0 ./ 2)
end

###########################################################################################################

struct ObsGreenVillainTime <: Obs
  obs_data::ObsData{Array{Complex{Float64},1}}
  green::Array{Complex{Float64},1}
  ObsGreenVillainTime(num_of_measure, Lt) = new(
    ObsDataArrayComplex("Green Villain Time", num_of_measure, ntuple(x -> Lt, 1)),
    zeros(ntuple(x -> Lt, 1)),
  )
end

function measure_obs!(obs_green::ObsGreenVillainTime, sim)
  fill!(obs_green.green, 0.0 + im * 0.0)
  for t = 1:sim.sim_params.Lt
    for j in eachindex(sim.lat.angle)
      j_site = Tuple(CartesianIndices(sim.lat.angle)[j])
      time_separation = t - 1
      new_time = mod1(j_site[3] + time_separation, sim.sim_params.Lt)
      obs_green.green[t] +=
        exp.(
          im * (
            sim.lat.angle[j_site[1], j_site[2], new_time] -
            sim.lat.angle[j_site...]
          )
        )
    end
  end
  take_measurement!(obs_green.obs_data, obs_green.green)
  return nothing
end

###########################################################################################################

#single particle Green's function for villain part
struct ObsGreenVillainSpace <: Obs
  obs_data::ObsData{Array{Complex{Float64},1}}
  green::Array{Complex{Float64},1}
  ObsGreenVillainSpace(num_of_measure, L) = new(
    ObsDataArrayComplex("Green Villain Space", num_of_measure, ntuple(x -> L, 1)),
    zeros(ntuple(x -> L, 1)),
  )
end

function measure_obs!(obs_green::ObsGreenVillainSpace, sim)
  fill!(obs_green.green, 0.0 + im * 0.0)
  for r in 1:sim.sim_params.L
    for j in eachindex(sim.lat.angle)
      separation = r - 1
      j_site = Tuple(CartesianIndices(sim.lat.angle)[j])
      new_x = mod1(j_site[1] + separation, sim.sim_params.L)
      obs_green.green[r] +=
        exp.(
          im *
          (sim.lat.angle[new_x, j_site[2], j_site[3]] - sim.lat.angle[j_site...])
        ) / 2
      new_y = mod1(j_site[2] + separation, sim.sim_params.L)
      obs_green.green[r] +=
        exp.(
          im *
          (sim.lat.angle[j_site[1], new_y, j_site[3]] - sim.lat.angle[j_site...])
        ) / 2
    end
  end
  take_measurement!(obs_green.obs_data, obs_green.green)
  return nothing
end

function cal_mag_villain_list(sim, time::Int64)
  total_mag = 0.0 + im * 0.0
  for i in 1:size(sim.lat.angle, 1)
    for j in 1:size(sim.lat.angle, 2)
      total_mag += exp.(im * sim.lat.angle[i, j, time])
    end
  end
  return total_mag
end

###########################################################################################################
