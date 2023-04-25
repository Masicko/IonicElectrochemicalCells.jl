############# This still produces equal charges L R ##########
function example_symmetry_charge()
    Rcell = AYALGBoltzmann()
    top_prms_dict = Dict(
        :AYSZ => 0.4799560546875, 
        :alpha => 1.0,
        #:alpha => 0.09976136948679738, 
    );
    update_parameters!(Rcell, top_prms_dict)
    Rcell.U = inival(Rcell)
	Rcell.Uold = copy(Rcell.U)
    stationary_update!(Rcell, Dict([]), tend=1e-5) 
    return Rcell
end



######################### EXAMPLES ################################


function get_my_top_AYALGBoltzmann()
    Rcell = AYALGBoltzmann()
    top_prms_dict = Dict(
        :AYSZ => 0.4799560546875, 
        :alpha => 1.0,
        #:alpha => 0.09976136948679738, 
    ); 
    update_parameters!(Rcell, top_prms_dict)
    Rcell.U = inival(Rcell)
	Rcell.Uold = copy(Rcell.U)
    stationary_update!(Rcell, Dict([]), tend=1e-5)   
    return Rcell
end



function get_my_top_AYALG1iBoltzmann(grid=AYALG1iBoltzmann().system.grid)
    top_prms_dict = Dict(
        #:DYSZ => 10e-9,
        :kA => 1.0e30,
        :AYSZ => 0.0,
        :alpha => 0.2, 
        :alphas => 0.5,
        :GA => 0.1*e0,
        :AYSZs => 0.0,
        :Ge => 0.0*e0,
        );
    Rcell = AYALG1iBoltzmann(grid)
	update_parameters!(Rcell, top_prms_dict)
    Rcell.U = inival(Rcell)
	Rcell.Uold = copy(Rcell.U)
    stationary_update!(Rcell, Dict([]), tend=1e-5)
    return Rcell
end

function get_my_top_AYALG1iBoltzmann_HALF(grid=AYALG1iBoltzmann_HALF().system.grid)
    top_prms_dict = Dict(
        #:DYSZ => 10e-9,
        :kA => 1.0e30,
        :AYSZ => 0.0,
        :alpha => 0.2, 
        :alphas => 0.5,
        :GA => 0.1*e0,
        :AYSZs => 0.0,
        :Ge => 0.0*e0,
        );
    Rcell = AYALG1iBoltzmann_HALF(grid)
	update_parameters!(Rcell, top_prms_dict)
    Rcell.U = inival(Rcell)
	Rcell.Uold = copy(Rcell.U)
    stationary_update!(Rcell, Dict([]), tend=1e-5)
    return Rcell
end

function get_my_top_prms_cell()
    #get_my_top_AYALGBoltzmann()
    get_my_top_AYALG1iBoltzmann()
end

L_R_charge(cell) = (
			get_charge(cell, "l")
)

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

function chargeshow(;bend = 1.0, bstep =0.05, tend=1e+0)
	Rcell = get_my_top_prms_cell()
    df = DataFrame(bias = [], U = [], charge = [])
	for bias ∈ collect(0.0 : bstep : bend)
		stationary_update!(Rcell, Dict(:bias => bias), tend=tend)
        push!(df, 
            (bias, copy(Rcell.U), L_R_charge(Rcell))
        )
	end 
    Rcell = get_my_top_prms_cell()
    for bias ∈ collect(-bstep : -bstep : -bend)
		stationary_update!(Rcell, Dict(:bias => bias), tend=tend)
        push!(df, 
            (bias, copy(Rcell.U), L_R_charge(Rcell))
        )
	end
    return sort!(df)
end

function eval_capacitance!(bdf)
    
    dbias = bdf[!, :bias][2] - bdf[!, :bias][1]

    capacitance_storage = [-Inf]
    for row in 2:size(bdf,1)-1
        Qs = bdf[:, :charge][row-1:row+1]
        cap= @. ((Qs[2] - Qs[1])/(dbias) + (Qs[3] - Qs[2])/(dbias))/2
        push!(capacitance_storage, cap)
    end
    push!(capacitance_storage, -Inf)
    bdf[!, :capacitance] = capacitance_storage
  
    return bdf
end

function test_capacitance(;bend = 1.0, bstep=0.05, tend=1e-0)
    bdf = chargeshow(;bend = bend +bstep, bstep =bstep, tend=tend)
    eval_capacitance!(bdf)
    return bdf
end

function plot_capacitance(p, bdf)
    plot!(p, bdf[!, :bias], bdf[!, :capacitance])
    return p
end




function chargeshow!(cell::AYALG1iBoltzmann; bend = 1.0, bstep=0.05, tend=1.0e-2)
    aux_df = DataFrame(bias = [], charge = [])
    function compute_charge(bias, cell)
        push!(aux_df, (bias, get_stored_charge(cell)))
        return
    end
    gdf = biasshow!(cell, bend=bend, bstep=bstep, tend=tend, callback=compute_charge)
    sort!(aux_df, [:bias])
    gdf[:, :charge] = aux_df[:, :charge]
    return gdf
end

function capacitance_measurement!(cell::AYALG1iBoltzmann;bend = 1.0, bstep=0.01, tend=1.0e-2)
    bdf = chargeshow!(cell, bend=bend+bstep, bstep=bstep, tend=tend)
    eval_capacitance!(bdf)
    delete!(bdf, [1, size(bdf)[1]])
    return bdf
end