using LinearAlgebra
using Statistics


# Obsevables
abstract type Obs end

mutable struct ObsData{T}
  name::String
  mean::T #direct sum
  square::T #sum of square for computing variance
  num_of_measure::Int64 #Total number of measurements
  c_num_of_measure::Int64 #Current number of measurements
end

########## Data type ##########################################################################################################################################################

#scalar
function ObsDataScalar(name::String, num_of_measure::Int64)
  ObsData{Float64}(name, 0.0, 0.0, 0, 0)
end

#fixed size real array, e.g. correlation function
function ObsDataArray(name::String, num_of_measure::Int64, dims::NTuple{D,Int64}) where {D}
  ObsData{Array{Float64,D}}(name, zeros(dims), zeros(dims), 0, 0)
end

#fixed size complex array, i.e. Fourier transform of stiffness
function ObsDataArrayComplex(
  name::String,
  num_of_measure::Int64,
  dims::NTuple{D,Int64},
) where {D}
  ObsData{Array{Complex{Float64},D}}(
    name,
    zeros(dims) .+ im * zeros(dims),
    zeros(dims) .+ im * zeros(dims),
    0,
    0,
  )
end

#Extensible real list
function ObsDataList(name::String, num_of_measure::Int64)
  ObsData{Array{Float64,1}}(name, zeros(num_of_measure), zeros(num_of_measure), 0, 0)
end

#Extensible complex list
function ObsDataComplex(name::String, num_of_measure::Int64)
  ObsData{Complex{Float64}}(name, 0.0 + im * 0.0, 0.0 + im * 0.0, 0, 0)
end

#Extensible list for quahtum hall conductivity
function ObsDataHall(name::String, num_of_measure::Int64)
  ObsData{Array{Float64,1}}(name, zeros(num_of_measure), zeros(num_of_measure), 0, 0)
end


########### Take Measurement ########################################################################################################################################################


#scalar
function take_measurement!(obs_data::ObsData{T}, data::T) where {T<:Number}
  obs_data.c_num_of_measure += 1
  obs_data.mean += data
  obs_data.square += data^2
  return nothing
end


#list, e.g. like green function
function take_measurement!(obs_data::ObsData{Array{T,D}}, data::Array{T,D}) where {T,D}
  obs_data.c_num_of_measure += 1
  for i in eachindex(data)
    obs_data.mean[i] += data[i]
    obs_data.square[i] += data[i]^2
  end
  return nothing
end

function clear_data!(measurement::T) where {T<:Obs}
  measurement.obs_data.mean = zero(measurement.obs_data.mean)
  measurement.obs_data.square = zero(measurement.obs_data.square)
  measurement.obs_data.c_num_of_measure = 0
  return nothing
end

########## Measurement methods ###################################################################################################################################################################

mutable struct MeasurementArray
  measurements::Array{Any,1}
  MeasurementArray(measurements::Array{Any,1}) = new(measurements)
end

function measure_all!(measurements::MeasurementArray, sim)
  measure_obs!.(measurements.measurements, Ref(sim))
  return nothing
end

function measure_correlator!(measurements::MeasurementArray, sim)
  measure_obs!.(measurements.measurements[1:3], Ref(sim))
  return nothing
end

function reset_measurements!(measurements::MeasurementArray)
  clear_data!.(measurements.measurements)
  return nothing
end

function get_measurements(m_array::MeasurementArray)
  """
  Return a dictionary of all observables
  """
  return Dict(
    m.obs_data.name => (m.obs_data.mean, m.obs_data.square) for
    m in m_array.measurements
  )
end
