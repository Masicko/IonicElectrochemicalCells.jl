using IonicElectrochemicalCells; iec = IonicElectrochemicalCells
using Test

@testset "IonicElectrochemicalCells.jl" begin
    # YHElectrolyte
    # halfelectrolytegrid = half_electrolyte_1D(1.0,4e-3)
    # cell = YHElectrolyte(halfelectrolytegrid)
    # biastest!(cell)
    # @test isapprox(IonicElectrochemicalCells.impedance_current(cell.system, cell.U), -0.3368816004135132)
    
    # # AYABoltzmann, AYALGBoltzmann
    # tgrid  = simplefullcell1D(1e-3, 1.0e-5; N=5)
    # testvals = Dict(
    #     :AYABoltzmann => -0.0, 
    #     :AYALGBoltzmann => -5.878539927672726e8,
    #     :AYALG1iBoltzmann => -1.0,
    # )
    # for cellcstr in [AYABoltzmann, AYALGBoltzmann]
    #     cell = cellcstr(tgrid)
    #     biastest!(cell)
    #     @test isapprox(get_stored_charge(cell), testvals[Symbol(cellcstr)])
    # end
    # cell = AYALG1iBoltzmann() # TODO use smaller grid > faster test BUT convergence...
    # biastest!(cell)
    # @test isapprox(get_stored_charge(cell),-0.37091832870777824)  # check the value in the DL 
    @test maximum(abs.(iec.test_half_cell())) < 1.0e-2
end
