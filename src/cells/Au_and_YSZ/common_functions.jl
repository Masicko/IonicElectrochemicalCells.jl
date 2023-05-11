# general functions for each cell
function update_parameters!(sys::VoronoiFVM.AbstractSystem, params_dict::Dict)
    set_parameters!(sys.physics.data, params_dict)
    boundary_dirichlet!(sys, ipsi, Γ_Aul, sys.physics.data.bias)
end

update_parameters!(cell::AbstractCell, params_dict::Dict) = update_parameters!(cell.system, params_dict)


function stationary_update!(cell::AbstractCell, params_dict; tend=1e-3)
    reload = false
    for par in [:alpha, :alphas, :SL, :SR]
        if (par ∈ keys(params_dict))
            if (getfield(cell.system.physics.data, par) != params_dict[par])
                reload = true
            end
        end
    end
    update_parameters!(cell, params_dict)
    cell.Uold .= cell.U
    if reload
        cell.U = inival(cell) # no re-allocations
    end
    tsol = VoronoiFVM.solve(
        cell.U,
        cell.system,
        [0.0, tend],
        control=control
    )
    cell.U .= tsol[end]
end

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
