mutable struct Lattice
    #angle and P variable of Villain model
    angle::Array{Float64,NDIMS}
    p_lattice::Array{Int64,NDIMS + 1}
    #current variable of the loop model
    current::Array{Int64,NDIMS + 1}
    #curl lattice is the lattice of the curl of the p variable
    curl_lattice::Array{Int64,NDIMS + 1}
    #fluctuating boundary
    fluctuate_lattice::Array{Float64,1}
    # head and tail
    head::NTuple{NDIMS,Int64}
    tail::NTuple{NDIMS,Int64}
    # winding number 
    winding::Array{Int64,1}
    energy::Float64

    function Lattice(sim_params::SimParams)
        L = sim_params.L
        T = sim_params.Lt
        angle_lattice = zeros(Float64, L * ones(Int64, NDIMS - 1)..., T)
        fluctuate_lattice = zeros(Float64, NDIMS)
        p_lattice = zeros(Int64, L * ones(Int64, NDIMS - 1)..., T, NDIMS)
        current_lattice = zeros(Int64, L * ones(Int64, NDIMS - 1)..., T, NDIMS)
        curl_lattice = zeros(Int64, L * ones(Int64, NDIMS - 1)..., T, NDIMS)
        new(
            angle_lattice,
            p_lattice,
            current_lattice,
            curl_lattice,
            fluctuate_lattice,
            tuple(ones(Int64, NDIMS)...),
            tuple(ones(Int64, NDIMS)...),
            zeros(Int64, NDIMS),
            0.0, #energy
        )
    end
end
