using IonicElectrochemicalCells; iec = IonicElectrochemicalCells
using Test


function AYA_Sl_biaswise_V_ISR_errors()
    STcell = AYA_Sl(iec.doublehalf_cell1D(
        iec.electrode_thickness, 
        iec.electrolyte_thickness, 
        iec.electrode_thickness, 
        dmin=1e-14))
    #STcell = AYA_Sl()
    stationary_update!(STcell,Dict([
				:AYSZ => 0.0,
                :AYSZs => 0.0,
				
				:GA => 0.1*e0,
				:Ge => 0.1*e0,
				:alpha => 0.001,
                :alphas => 0.001,

                :SL => 0.5,
                :SR => 0.7,
				
                :DYSZ => 1.0e-8,
				:kA => 1.0e25,

				#:kA => 0.0,
				#:boundary_charge_fac => 0.0
			]))

    return maximum(abs.(iec.biaswise_V_ISR_errors(STcell))) < 1.0e-2
end

@testset "IonicElectrochemicalCells.jl" begin
    @test AYA_Sl_biaswise_V_ISR_errors()
end