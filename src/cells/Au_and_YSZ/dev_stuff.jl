######################### EXAMPLES ################################


function U_bias_symmetry_test!(bdf)
    symm_column = []
    Us = bdf[!, :U]
    for i in 1:Int((length(Us)-1)/2)
        
        push!(symm_column,
            maximum(abs.(Us[i][1,:] .- reverse(Us[end+1-i][1,:] .+ bdf.bias[i])))
        )
    end
    #the bias 0.0K
    push!(symm_column, 0.0)
    append!(symm_column, reverse(symm_column[1:end-1]))
    bdf[!, :U_bsymmetry] = symm_column
    
    return bdf 
end

function testing_two_cell(;my_bias=0.1, tend=1e+0)
	Rcell = get_my_top_prms_cell()
	Lcell = get_my_top_prms_cell()
	for bias ∈ collect(0.0 : 0.1 : my_bias)
		stationary_update!(Rcell, Dict(:bias => bias), tend=tend)
		stationary_update!(Lcell, Dict(:bias => -bias), tend=tend)
	end
	L_R_charge(cell) = (
			IonicElectrochemicalCells.get_charge(cell, "l") 
			+
			IonicElectrochemicalCells.get_charge(cell, "r")
		)/2
	@show L_R_charge(Rcell)
	@show L_R_charge(Lcell)

    @show IonicElectrochemicalCells.get_charge(Rcell, "l")
    @show IonicElectrochemicalCells.get_charge(Rcell, "r") 

    @show IonicElectrochemicalCells.get_charge(Lcell, "l")
    @show IonicElectrochemicalCells.get_charge(Lcell, "r") 
    return Lcell, Rcell
end

function test_capacitance(;bend = 1.0, bstep=0.05, tend=1e-0)
    bdf = chargeshow(;bend = bend +bstep, bstep =bstep, tend=tend)
    eval_capacitance!(bdf)
    return bdf
end


function get_comparable_quantities(STcell)
	if typeof(STcell) == IEC_const.AYA_Sl
        iec = eval(Symbol(:IEC_const))
    elseif typeof(STcell) == IonicElectrochemicalCells.AYA_Sl
        iec = eval(Symbol(:IonicElectrochemicalCells))
    else
        println("ERROR: wrong type incomiiiing")
        return
    end
    data = STcell.system.physics.data
	X = STcell.system.grid.components[ExtendableGrids.Coordinates]
	BFNodes = STcell.system.grid.components[ExtendableGrids.BFaceNodes]
	testing_id = b1_id = BFNodes[3]
	U = STcell.U
	U_ISR = U[:, b1_id]

    b2_id = BFNodes[4]
    X_list = X[b1_id : b2_id]
    center = Int(round(length(X_list)/2))
    phi_C = U[1, b1_id+center]
	
	yV = U[2, testing_id]
	phi_ISR = U[1, testing_id]
	V_ISR = phi_ISR - phi_C
	dphi = U[1, testing_id+1] - U[1, testing_id]
	dx = X[testing_id+1]-X[testing_id]
	println("...")
    @show V_ISR
	
	println("\n <<<< Comparable quantities to analytical solution: >>>>")
	println(" --- YSZ side:")
	@show yV dphi/dx

	phi_L = U[1, 1]
	U_Au = U[1, testing_id] - phi_L
	ye_eq = iec.BoltzmannAu_ne(data, U_Au)/iec.nLAu
	dphi_Au = U[1, testing_id] - U[1, testing_id-1]
	dx_Au = X[testing_id] - X[testing_id-1]
	println(" --- Au side:")
	@show ye_eq dphi_Au/dx_Au

	
	nVs_eq = iec.nVmaxs(data.alphas, data.SL ) * U[3, testing_id]
	nes_eq = iec.ISR_electrondensity(U_ISR, iec.dummy_bnode(iec.:Γ_YSZl), data)
	yVs_eq = U[3, testing_id]
	yes_eq = (1/iec.nLAus(data.SL))*nes_eq

	nFs = iec.ISR_chargedensity(nes_eq, nVs_eq, data.SL)
	nFs_sim(nes, nVs) = iec.ISR_chargedensity(nes, nVs, data.SL)
	println(" --- ISR:")
	@show yVs_eq yes_eq nFs
end

# using IEC_const

# function both_get_cells(prms_dict)
#     Ccell = IEC_const.AYA_Sl()
#     IEC_const.update_parameters!(Ccell, prms_dict)
#     Ccell.U = IEC_const.inival(Ccell)
#     IEC_const.stationary_update!(Ccell, Dict(:bias => 0.0))

#     STcell = AYA_Sl()
#     update_parameters!(STcell, prms_dict)
#     STcell.U = inival(STcell)
#     stationary_update!(STcell, Dict(:bias => 0.0))
#     return STcell, Ccell
# end

# function both_stationary_update!(STcell, Ccell, prms_dict)
#     IEC_const.update_parameters!(Ccell, prms_dict)
#     IEC_const.stationary_update!(Ccell, prms_dict, tend=1.0e2)
#     update_parameters!(STcell, prms_dict)
#     stationary_update!(STcell, prms_dict, tend=1.0e2)
# end

# function test_consistency()
#     std_dict = Dict(
# 			:T => 800.0,

# 			:AYSZ => 0.0,
# 			#:alpha => 8.0e-2,
# 			:alpha => 8.0e-5,
			
# 			:AYSZs => 0.0,
# 			#:alphas => 0.025,
# 			:alphas => 0.0025,
# 			:GA => -0.1*e0,
# 			:Ge => -0.0*e0,
# 			:SL => 1.0,
# 			:SR => 1.0,
# 			:DYSZ => 1.0e-9,
# 			#:DYSZ => 1.0e-7,
		
			
# 			#:boundary_charge_fac => 0.0, 
# 			#:kA => 1.0e0
# 	)
#     STcell, Ccell = both_get_cells(std_dict);
#     stationary_update!(STcell, Dict(:SL => 0.5, :SR => 0.5))
#     both_stationary_update!(STcell, Ccell, Dict(:bias => 0.2))
    
#     get_comparable_quantities(STcell)
#     get_comparable_quantities(Ccell)

#     show_partial_charges(get_partial_charges(STcell))
#     show_partial_charges(IEC_const.get_partial_charges(Ccell))
#     return STcell, Ccell   
# end
