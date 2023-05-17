using IonicElectrochemicalCells; iec = IonicElectrochemicalCells
using Test

@testset "IonicElectrochemicalCells.jl" begin
    @test maximum(abs.(iec.test_AYA_Sl_cell())) < 3.0e-3
end
