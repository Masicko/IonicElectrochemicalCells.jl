bulk_charge(cell::AbstractCell) = -1 .*VoronoiFVM.integrate(cell.system, cell.system.physics.reaction, cell.U) # nspec x nregion
boundary_charge(cell::AbstractCell) = -1 .*VoronoiFVM.integrate(cell.system, cell.system.physics.breaction, cell.U, boundary=true)

function nofunc(bias, cell)
end

function biasshow!(cell::AbstractCell; bend=1.0, bstep=0.1, tend=1e-3, callback=nofunc)
    cell.U = inival(cell)
    cell.Uold = copy(cell.U)
    results = DataFrame(bias=Float64[], U=Any[])
    for bias ∈ collect(0.0:bstep:bend)
        stationary_update!(cell, Dict(:bias => bias), tend=tend)
        callback(bias, cell)
        push!(results, Dict(:bias => bias, :U => copy(cell.U)))
    end

    cell.U = inival(cell)
    cell.Uold = copy(cell.U)
    stationary_update!(cell, Dict(:bias => 0.0), tend=tend)
    for bias ∈ collect(-bstep:-bstep:-bend)
        stationary_update!(cell, Dict(:bias => bias), tend=tend)
        callback(bias, cell)
        push!(results, Dict(:bias => bias, :U => copy(cell.U)))
    end
    return sort!(results, [:bias])
end


function stateplot(cell::AbstractCell, U; Plotter=nothing, Xscale=false,xlim=nothing, ylim=nothing)
	if Xscale == true
		gr = cell.system.grid.components[Coordinates][:]
	else
		gr = 1:length(U[1,:])
	end
	if xlim!=nothing
        if ylim!=nothing
            p = Plotter.plot(title = "Coverage", xlimits=xlim, limits=ylim)
        else
		    p = Plotter.plot(title = "Coverage", xlimits=xlim)
        end
	else
        if ylim!=nothing
            p = Plotter.plot(title = "Coverage", limits=ylim)
        else
		    p = Plotter.plot(title = "Coverage")
        end
	end
    for i in 1:size(U)[1]
        Plotter.plot!(p, gr, U[i,:], linewidth = 2, markertype=:circle, markersize=5)
    end
	return p
end

function stateplot(cell::Union{AYALG1iBoltzmann_HALF,AYALG1iBoltzmann, AYA_Sl} , U; 
    Plotter=nothing, Xscale=false,xlim=nothing, ylim=nothing,
    title = "Coverage"
    )
	phi_L = U[1,1]
    data = cell.system.physics.data
    phi_to_ye(V) = BoltzmannAu_ne(data, (V)) / nLAu

    ye_data_L = phi_to_ye.(U[1,:] .- phi_L)
    boundary_id = cell.system.grid.components[ExtendableGrids.BFaceNodes][3]

    if Xscale == true
		gr = cell.system.grid.components[Coordinates][:]
	else
		gr = 1:length(U[1,:])
	end
    if xlim!=nothing
        if ylim!=nothing
            p = Plotter.plot(title = title, xlimits=xlim, limits=ylim)
        else
		    p = Plotter.plot(title = title, xlimits=xlim)
        end
	else
        if ylim!=nothing
            p = Plotter.plot(title = title, limits=ylim)
        else
		    p = Plotter.plot(title = title)
        end 
	end

    Plotter.plot!(p, gr, U[1,:], linewidth = 2, label="phi")
    Plotter.plot!(p, gr, U[2,:], linewidth = 2, label="yV")
    Plotter.scatter!(p, [gr[boundary_id]], [U[3,boundary_id]], label="yVs_L")
    Plotter.plot!(p, gr[1:boundary_id], ye_data_L[1:boundary_id], label="ye_L")

    if typeof(cell) in [AYALG1iBoltzmann, AYA_Sl]
        ye_data_R = phi_to_ye.(U[1,:])
        boundary2_id = cell.system.grid.components[ExtendableGrids.BFaceNodes][4]
        Plotter.plot!(p, gr[boundary2_id:end], ye_data_R[boundary2_id:end], label="ye_R")
        Plotter.scatter!(p, [gr[boundary2_id]], [U[3,boundary2_id]], label="yVs_R")
    end
	return p
end

##### charge evaluation

function chargeshow!(cell::Union{AYALG1iBoltzmann, AYA_Sl}; bend = 1.0, bstep=0.05, tend=1.0e-2, testing=false)
    if testing
        aux_df = DataFrame(bias = [], charge = [], V_ISR_error = [])
    else
        aux_df = DataFrame(bias = [], charge = [])
    end
    function compute_charge(bias, cell)
        if testing
            if !test_partial_charges(cell) 
                println("bias: ",bias,"  -> TEST_PARTIAL_CHARGES  ... FAILED!")
            end
            push!(aux_df, (bias, get_stored_charge(cell), get_V_ISR_error(cell)))
        else
            push!(aux_df, (bias, get_stored_charge(cell)))
        end
        return
    end
    gdf = biasshow!(cell, bend=bend, bstep=bstep, tend=tend, callback=compute_charge)
    sort!(aux_df, [:bias])
    gdf[:, :charge] = aux_df[:, :charge]
    if testing
        gdf[:, :V_ISR_error] = aux_df[:, :V_ISR_error]
    end
    return gdf
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

function capacitance_measurement!(cell::Union{AYALG1iBoltzmann, AYA_Sl};bend = 1.0, bstep=0.01, tend=1.0e-2, testing=false)
    bdf = chargeshow!(cell, bend=bend+bstep, bstep=bstep, tend=tend, testing=testing)
    eval_capacitance!(bdf)
    delete!(bdf, [1, size(bdf)[1]])
    return bdf
end

function plot_capacitance(p, bdf)
    plot!(p, bdf[!, :bias], bdf[!, :capacitance])
    return p
end

function show_partial_charges(pch)
    println("  <<<<   ===========   partial charges  ==========  >>>>")
    println("Left Au:   ", pch[:nF_Au_L] + pch[:bQ_Au_L])
    println("Left YSZ:  ", pch[:nF_YSZ_L] + pch[:bQ_YSZ_L])
    println("Right YSZ: ", pch[:nF_YSZ_R] + pch[:bQ_YSZ_R])
    println("Right Au:  ", pch[:nF_Au_R] + pch[:bQ_Au_R])
    println("")
    println("Left DL:   ", pch[:nF_Au_L] + pch[:bQ_Au_L] + pch[:bQ_YSZ_L] + pch[:nF_YSZ_L])
    println("Right DL:  ", pch[:nF_Au_R] + pch[:bQ_Au_R] + pch[:bQ_YSZ_R] + pch[:nF_YSZ_R])
    println("")
    println("total YSZ: ", pch[:bQ_YSZ_L] + pch[:nF_YSZ_L] + pch[:nF_YSZ_R] + pch[:bQ_YSZ_R])
    return
end