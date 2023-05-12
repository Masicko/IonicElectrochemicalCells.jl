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