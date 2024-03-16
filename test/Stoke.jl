#check if Stoke theorem is obeyed, i.e. curl_P = etagral of loop in P variable around curl_P vector
function LoopIntegral(c_site::NTuple{Defs.NDIMS,Int64}, dir1::Int64, sim::MCMC.Sim)
  sum = 0
  dir2 = Defs.cyclic[dir1+1]
  dir3 = Defs.cyclic[dir1+2]

  loop_site = MCMC.hop(c_site, dir1, 1, size(sim.lat.angle))

  sum += sim.lat.p_lattice[loop_site..., dir2]
  sum -= sim.lat.p_lattice[loop_site..., dir3]

  jump1 = MCMC.hop(loop_site, dir2, 1, size(sim.lat.angle))
  sum += sim.lat.p_lattice[jump1..., dir3]

  jump2 = MCMC.hop(loop_site, dir3, 1, size(sim.lat.angle))
  sum -= sim.lat.p_lattice[jump2..., dir2]
  return sum
end

@testset "Stoke" begin
  for i in eachindex(sim.lat.angle)
    c_site = Tuple(CartesianIndices(sim.lat.angle)[i])
    for dir = 1:Defs.NDIMS
      #now integrate/sum along the loops
      @test LoopIntegral(c_site, dir, sim) == sim.lat.curl_lattice[c_site..., dir]
    end
  end
end
