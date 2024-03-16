@testset "Angle" begin
  @testset "Angle lattice" begin
    for site in eachindex(sim.lat.angle)
      site = Tuple(CartesianIndices(sim.lat.angle)[site])
      @test -pi < sim.lat.angle[site...] < pi
    end
  end

  @testset "Fluctuate Angle" begin
    for dir in 1:Defs.NDIMS
      @test -pi < sim.lat.fluctuate_lattice[dir] < pi
    end
  end
end
