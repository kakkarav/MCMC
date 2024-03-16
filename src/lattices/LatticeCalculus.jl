using Base.Cartesian

@generated function hop(
  index::NTuple{NDIMS,Int64},
  dir::Int64,
  lr::Int64,
  dims::NTuple{NDIMS,Int64},
)
  quote
    @nif 2 lr_d -> (lr == lr_d) lr_d -> begin
      @nif $NDIMS d -> (d == dir) d -> begin
        @ntuple $NDIMS d_t -> begin
          if d_t == d
            if (lr_d == 1)
              (index[d_t] == dims[d_t] ? 1 : index[d_t] + 1)
            elseif (lr_d == 2)
              (index[d_t] == 1 ? dims[d_t] : index[d_t] - 1)
            else
              print(" lr is out of bound")
            end
          else
            index[d_t]
          end
        end
      end
    end
  end
end

function curl(lattice::Lattice, location::NTuple{NDIMS,Int64}, dir::Int64)
  villain_site = hop(location, dir, 1, size(lattice.angle))
  n = cyclic[dir+1]
  nn = cyclic[dir+2]
  return derivative(lattice, villain_site, n, nn) -
         derivative(lattice, villain_site, nn, n)
end

function derivative(
  lattice::Lattice,
  location::NTuple{NDIMS,Int64},
  dir::Int64,
  value::Int64,
)
  nn = hop(location, dir, 1, size(lattice.angle))
  return lattice.p_lattice[nn..., value] - lattice.p_lattice[location..., value]
end
