using TSMC: Defs, SimParameters, MCMC

@testset "Curl" begin
    for site in eachindex(sim.lat.angle)
        site = Tuple(CartesianIndices(sim.lat.angle)[site])
        for dir in 1:Defs.NDIMS
            @test MCMC.curl(sim.lat, site, dir) == sim.lat.curl_lattice[site..., dir]
        end
    end
end
