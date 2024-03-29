module MCMC

# Include parameter definitions and simulation parameters
include("parameters/defs.jl")
include("parameters/simParameters.jl")

using .Defs, .SimParameters, Statistics, Random

# Include lattice data and lattice calculus
include("lattices/latticedata.jl")
include("lattices/disorderdata.jl")
include("lattices/latticeCalculus.jl")

# Include measurement and update modules
include("measurements/measurements.jl")
include("measurements/observablesScalar.jl")
include("measurements/observablesList.jl")
include("updates/updateAngle.jl")
include("updates/updateP.jl")
include("updates/updateJ.jl")
include("updates/updateFluctuate.jl")
include("updates/updateBinding.jl")
include("updates/updateCustom.jl")



struct Sim
  lat::Lattice
  lattice_map::CartesianIndices{NDIMS}
  worm_data::WormMoveData
  sim_params::SimParams
  measurements::MeasurementArray
  rng::Random.Xoshiro
  fourier::Array{Array{Complex{Float64},1},1}
  disorder::Disorder
end

function Sim(sim_params::SimParams)
  rng = Random.Xoshiro(rand(UInt64))
  fourier = [
    [exp(-im * 2 * pi * (n - 1) * (m - 1) / sim_params.Lt) for n in 1:sim_params.Lt] for m in 1:sim_params.Lt
  ]
  lattice = Lattice(sim_params)
  lattice_map = CartesianIndices(lattice.angle)
  Sim(
    lattice,
    lattice_map,
    WormMoveData(sim_params, rng),
    sim_params,
    MeasurementArray(
      Any[
        ObsGreenLoopTime(sim_params.num_of_measure, sim_params.Lt),
        ObsGreenLoopSpace(sim_params.num_of_measure, sim_params.L),
        ObsZRatio(sim_params.num_of_measure),
        ObsCompressLoop(sim_params.num_of_measure),
        ObsStiffnessLoopAll(sim_params.num_of_measure, sim_params.Lt),
        ObsEnergy(sim_params.num_of_measure),
        ObsHallConductivity(sim_params.num_of_measure, sim_params.Lt),
        ObsCurrentTimeLoop(sim_params.num_of_measure),
        ObsCurrentSpaceLoop(sim_params.num_of_measure),
        ObsCurrentLoopFT(sim_params.num_of_measure, sim_params.Lt),
        ObsGreenVillainTime(sim_params.num_of_measure, sim_params.Lt),
        ObsGreenVillainSpace(sim_params.num_of_measure, sim_params.L),
        ObsCompressVillain(sim_params.num_of_measure),
        ObsCurrentSpaceVillain(sim_params.num_of_measure),
        ObsCurrentTimeVillain(sim_params.num_of_measure),
        ObsCurrentVillainFT(sim_params.num_of_measure, sim_params.Lt),
        ObsStiffnessVillain(sim_params.num_of_measure),
        ObsStiffnessVillainFT(sim_params.num_of_measure, sim_params.Lt),
        ObsMagnitization(sim_params.num_of_measure),
      ],
    ),
    rng,
    fourier,
    Disorder(sim_params),
  )
end

function move_worm_head!(sim::Sim)
  """Implement the classical worm algorithm"""
  closed = (sim.lat.head == sim.lat.tail)
  if closed
    jump_or_shift = rand(sim.rng)
    if (0.5 > jump_or_shift)
      return jump_worm!(sim.worm_data, sim.lat)
    end
  end
  shift_worm!(sim.worm_data, sim.lat, sim.sim_params, sim.disorder)
end


function update_worm!(sim::Sim)
  for _ in 1:(3*sim.sim_params.L^3)
    move_worm_head!(sim)
  end
  return nothing
end

function update_villain!(sim::Sim)
  @inbounds for site in eachindex(sim.lat.angle)
    location = Tuple(sim.lattice_map[site])
    metro_angle!(sim.lat, sim.sim_params, sim.disorder, sim.rng, location)
  end
  @inbounds for dir in 1:NDIMS, site in eachindex(sim.lat.angle)
    location = Tuple(sim.lattice_map[site])
    metro_p!(sim.lat, sim.sim_params, sim.disorder, sim.rng, location, dir)
  end
  return nothing
end

function update_binding!(sim::Sim)
  #we update the villain part by sweeping through all bonds
  for dir in 1:NDIMS, site in eachindex(sim.lat.angle)
    location = Tuple(sim.lattice_map[site])
    metro_binding!(sim.lat, sim.sim_params, sim.disorder, sim.rng, location, dir)
  end
  return nothing
end

function update_all!(sim::Sim)
  update_villain!(sim)
  update_binding!(sim)
  update_worm!(sim)
  return nothing
end

function thermalize!(sim::Sim)
  for _ in 1:sim.sim_params.num_of_thermal
    update_all!(sim)
  end
  return nothing
end

function run!(sim::Sim)
  counter = 0
  while counter < sim.sim_params.num_of_measure
    for _ in 1:sim.sim_params.num_of_sweeps
      update_all!(sim)
    end
    measure_correlator!(sim.measurements, sim)
    if (sim.lat.head == sim.lat.tail)
      measure_all!(sim.measurements, sim)
      counter += 1
    end
  end
end

include("Helpers.jl")

end
