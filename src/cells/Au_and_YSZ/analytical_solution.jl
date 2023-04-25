function HALF_test_analytical(data, V_tot, V_ISR; verbose=true)
    save_sqrt(x) = ( x >= 0 ? sqrt(x) : ( x > -eps() ? 0.0 : sqrt(x) ))

    T = data.T

    # YSZ bulk
    yB = yV_en(data.alpha)
    X(V) = yB/(1-yB)*exp(-(zV*e0)/(kB*T)*V)

    yV_eq(V) = X(V)/(X(V) + 1)
    dphi_dx_YSZ(V) = -sign(V)*sqrt((2*e0*nLYSZ/epsYSZ))*save_sqrt(
        (
            -(mC*zC - mA*zV)*V
            +
            (kB*T)/(zV*e0)*
            (mA*zV - (1-data.alpha)*mC*zC)*
            log((1-yB)*(X(V) + 1))
        )
    )
    #cap_YSZ(V) = YSZ_charge_density(nVmax(data.alpha) * yV_eq(V)) / dphi_dx_YSZ(V)

    # Au bulk
    ye_eq(V) = exp(-(ze*e0)/(kB*T)*V)
    dphi_dx_Au(V) = sign(V)*sqrt((2*e0*nLAu/epsAu))*save_sqrt(
        -V
        +
        (kB*T)/(ze*e0)*(1 - ye_eq(V))
    )
    #cap_Au(V) = Au_charge_density(nLAu * ye_eq(V)) / dphi_dx_Au(V)
    
    # surface
    S = data.SL
    KV = exp(-data.GA / data.T / kB )
    nVs_eq(V, S)  = nVmaxs(data.alphas, S)* (KV*yV_eq(V))/(yV_eq(V)*(KV -1) + 1)
    nes_eq(V, S) = nLAus(S) * exp((-data.Ge - ze*e0*V)/(kB*T))
    #
    yVs_eq(V, S) = nVs_eq(V, S)/nVmaxs(data.alphas, S)
    yes_eq(V, S) = nes_eq(V, S)/nLAus(S)

    nFs(nes, nVs, S) = data.boundary_charge_fac*ISR_chargedensity(nes, nVs, S)
    
    poisson_eq(dphi_dx_Au, dphi_dx_YSZ, nVs, nes) = 
        -epsYSZ*dphi_dx_YSZ + epsAu*dphi_dx_Au - nFs(nes, nVs, S)

    U_Au(V_ISR) = V_ISR-V_tot
    U_YSZ(V_ISR)  = V_ISR
    
    poisson_eq(V_ISR) = poisson_eq(
        dphi_dx_Au(U_Au(V_ISR)),
        dphi_dx_YSZ(U_YSZ(V_ISR)),
        nVs_eq(U_YSZ(V_ISR), S),
        nes_eq(U_Au(V_ISR), S)
        )

    to_optimize(x) = poisson_eq(x[1])
    
    if verbose 
        println()
        println(" --- YSZ side:")
        @show yV_eq(U_YSZ(V_ISR)) dphi_dx_YSZ(U_YSZ(V_ISR))
        println(" --- Au side:")
        @show ye_eq(U_Au(V_ISR)) dphi_dx_Au(U_Au(V_ISR))
        println(" --- ISR:")
        @show yVs_eq(U_YSZ(V_ISR), S) yes_eq(U_Au(V_ISR), S)
    end
    return poisson_eq
end

# for the given voltage on a halfcell, it returns the semi-analytical 
# (solving numerically one algebraic equation) solution for potential on ISR
function HALF_get_analytic_V_ISR(data, V_tot)
    g = HALF_test_analytical(data, V_tot, 666, verbose=false)
    to_optimize(x) = g(x[1])^2
    Optim.optimize(to_optimize, -10, 10, GoldenSection()).minimizer
end

function test_half_cell()
    STcell = get_my_top_AYALG1iBoltzmann_HALF(
		 	half_cell1D(
		 		electrode_thickness, electrolyte_thickness, 	
				electrode_thickness, dmin=1.0e-18)
    )
    stationary_update!(STcell,Dict([
				#:DYSZ=>1.0e-9,
				
				:GA => 0.2*e0,
				:Ge => -0.2*e0,
				:AYSZ => 0.0,
				:alpha => 0.4,
                :alphas => 0.3,
				:kA => 1.0*1.0e25,
				
				#:kA => 0.0,
				#:boundary_charge_fac => 0.0
			]))

    function perform_it(bias, STcell)
        data = STcell.system.physics.data
        testing_id = STcell.system.grid.components[ExtendableGrids.BFaceNodes][3]
        phi_ISR = STcell.U[1, testing_id]
        phi_L = STcell.U[1, 1]
        phi_R = STcell.U[1, end-2]
        V_tot = phi_L - phi_R
        V_ISR = phi_ISR - phi_R
        #@show phi_R
        rel_phi_err = (HALF_get_analytic_V_ISR(data, V_tot) - V_ISR)/V_ISR
        #@show bias, rel_phi_err
        push!(rel_err_storage, rel_phi_err)
    end

    rel_err_storage = []
    biasshow!(STcell, bend=0.2, bstep=0.01, tend=1e-2, callback=perform_it)
    return rel_err_storage
end