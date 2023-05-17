
save_sqrt(x) = ( x >= 0 ? sqrt(x) : ( x > -eps() ? 0.0 : sqrt(x) ))

function HALF_test_analytical_AYA_Sl(data, V_tot, V_ISR; verbose=true, return_nFs=false)

    save_sqrt(x) = ( x >= 0 ? sqrt(x) : ( x > -eps() ? 0.0 : sqrt(x) ))

    #=
    ### YSZ bulk ###
    =#
    x_frac = 0.13
    epsYSZ = 28.0 * ε0          # https://link.springer.com/article/10.1007/BF01538778
    
    zCf(x) = 4*(1.0-x)/(1.0+x) + 2*3*(x/(1.0+x))
    zC = zCf(x_frac)
    mC = 4
    mA = 8
    zV = 2.0
    zi = -zV

    nL_YSZ_const = 7.46268656716418e27 # https://www.msesupplies.com/products/ysz-single-crystal?variant=31250580111418
    nL_YSZ(rS) = rS * nL_YSZ_const # maybe there should be some rS^2/3 or so.. but it does not matter now
    
    nF_YSZ(nV, rS) = e0*(nL_YSZ(rS)*(zC*mC - 2*mA) + zV*nV)
    nV_en(rS) = nL_YSZ(rS)*(mA - zC/zV*mC)
    nV_max(alpha, rS) = (1.0 - alpha)*nV_en(rS) + alpha*mA*nL_YSZ(rS)
    #yV(alpha, rS) = nV/nV_max(alpha, rS)  => 
    nV(yV, alpha, rS) = yV * nV_max(alpha, rS)
    #yV_en(rS, alpha) = nV_en(rS) / nV_max(alpha, rS)  # in fact, yV_en does NOT depend on rS
    yV_en(alpha) = (-mC * zC/zV + mA)/(-zC*mC / zV * (1 - alpha) + mA)

    nF_YSZ(yV, alpha, rS) = nF_YSZ(nV(yV, alpha, rS), rS)

    #=
    ### Au bulk ###
    =#

    a0  = 5.29177210903e-11       # m   Au lattice constant
    epsAu = 6.9 * ε0 # Separation of the contribution of free and bound electrons into real and imaginary parts of the dielectric constant of gold.     Shklyarevskii, I. N.; Pakhmov, P. L. USSR. Optika i Spektroskopiya  (1973),  34(1),  163-6.
    ze = -1.0
    nL_Au_const = 1 / (4 * pi / 3 * (3.01 * a0)^3) # RM Martin, Electronic structure (eq 5.1, )
    nL_Au(rS) = rS*nL_Au_const

    nF_Au(ne, rS) = e0 * (nL_Au(rS) * 1.0 + ze * ne)
    ne_en(rS) = nL_Au(rS)
    ne(ye, rS) = ye*nL_Au(rS)
    ye_en = 1.0

    #nF_Au(ye, rS) = nF_Au(ne(ye, rS), rS)
    # analytical solution to Au
    ye_Boltzmann(data, V) = exp(-ze * e0 / kB / data.T * V)
    ne_Boltzmann(rS, data, V) = nL_Au(rS) * ye_Boltzmann(data, V)
    
    #=
    ### ISR ###
    =#
    nLs_YSZ(rS::Float64) = (nL_YSZ(rS))^(2 / 3.0)  # [# YSZ faces/m^2]
    nLs_Au(rS::Float64) = (nL_Au(rS))^(2 / 3.0)    # [# Au faces/m^2]

    # YSZ part 
    mCs = 2               # [# cation sites/YSZ face]
    mAs = 4               # [# anion sites/YSZ face]

    nFs_YSZ(nVs, rS) = e0*(nLs_YSZ(rS)*(zC*mCs - 2*mAs) + zV*nVs)
    nVs_en(rS) = nLs_YSZ(rS)*(mAs - zC/zV*mCs)
    nVs_max(alphas, rS) = (1.0 - alphas)*nVs_en(rS) + alphas*mAs*nLs_YSZ(rS)
    #yVs(alphas, rS) = nVs/nVs_max(alphas, rS)  => 
    nVs(yVs, alphas, rS) = yVs * nVs_max(alphas, rS)
    yVs_en(alphas) = (-mCs * zC/zV + mAs)/(-mCs * zC/zV * (1 - alphas) + mAs)

    nFs_YSZ(yVs, alphas, rS) = nFs_YSZ(nVs(yVs, alphas, rS), rS)

    # Au part
    nFs_Au(nes, rS) = e0 * (nLs_Au(rS) * 1.0 + ze * nes)
    nes_en(rS) = nLs_Au(rS)
    nes(yes, rS) = yes*nLs_Au(rS)
    yes_en = 1.0

    ISR_arearatio(bnode,data) = (bnode.region == Γ_YSZl ? 
                                    data.SL : 
                                    (bnode.region == Γ_YSZr ? 
                                        data.SR : 
                                        println("ERROR: bnode.region = $(bnode.region)")
                                    )
                                )
    reduced_voltage(u, bnode, data) = u[ipsi] - (bnode.region == Γ_YSZl ? data.bias : 0.0)
    yes_Boltzmann(data, V) = exp(-data.Ge/ kB/ data.T) * ye_Boltzmann(data, V)  
    nes_Boltzmann(rS::Float64, data, V) = nLs_Au(rS)*yes_Boltzmann(data, V)
    nes_Boltzmann(u::Array, bnode, data) = nes_Boltzmann(ISR_arearatio(bnode, data), data, reduced_voltage(u,bnode,data))
    
    # common part
    nFs_tot(nes, nVs, rS) = nFs_YSZ(nVs, rS) + nFs_Au(nes, rS)
    nFs_tot(yes, yVs, alphas,  rS) = nFs_tot(nVs(yVs, alphas, rS), nes(yes, rS), rS)  




    ################################################################################################# 
    ################################### START Analytical solution ###################################
    ################################################################################################# 
    

    T = data.T
    rS = data.SL

    # YSZ bulk
    yB = yV_en(data.alpha)
    X(V) = yB/(1-yB)*exp(-(zV*e0)/(kB*T)*V)

    yV_eq(V) = X(V)/(X(V) + 1)
    dphi_dx_YSZ(V) = -sign(V)*sqrt((2*e0*nL_YSZ(rS)/epsYSZ))*save_sqrt(
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
    dphi_dx_Au(V) = sign(V)*sqrt((2*e0*nL_Au(rS)/epsAu))*save_sqrt(
        -V
        +
        (kB*T)/(ze*e0)*(1 - ye_eq(V))
    )
    #cap_Au(V) = Au_charge_density(nLAu * ye_eq(V)) / dphi_dx_Au(V)
    
    # surface
    KV = exp(-data.GA / data.T / kB )
    yVs_eq(V)  = (KV*yV_eq(V))/(yV_eq(V)*(KV -1) + 1)
    yes_eq(V) =  exp((-data.Ge - ze*e0*V)/(kB*T))
    #
    nVs_eq(V, rS)  = nVs_max(data.alphas, rS)* yVs_eq(V)
    nes_eq(V, rS) = nLs_Au(rS) * yes_eq(V)

    nFs(nes, nVs, rS) = data.boundary_charge_fac*nFs_tot(nes, nVs, rS)
    
    poisson_eq(dphi_dx_Au, dphi_dx_YSZ, nVs, nes) = 
        -epsYSZ*dphi_dx_YSZ + epsAu*dphi_dx_Au - nFs(nes, nVs, rS)

    U_Au(V_ISR) = V_ISR-V_tot
    U_YSZ(V_ISR)  = V_ISR
    
    poisson_eq(V_ISR) = poisson_eq(
        dphi_dx_Au(U_Au(V_ISR)),
        dphi_dx_YSZ(U_YSZ(V_ISR)),
        nVs_eq(U_YSZ(V_ISR), rS),
        nes_eq(U_Au(V_ISR), rS)
        )
    
    if verbose 
        println()
        println(" --- YSZ side:")
        @show yV_eq(U_YSZ(V_ISR)) dphi_dx_YSZ(U_YSZ(V_ISR))
        println(" --- Au side:")
        @show ye_eq(U_Au(V_ISR)) dphi_dx_Au(U_Au(V_ISR))
        println(" --- ISR:")
        @show yVs_eq(U_YSZ(V_ISR)) yes_eq(U_Au(V_ISR))
        @show nFs(  nes_eq(U_Au(V_ISR), rS),        nVs_eq(U_YSZ(V_ISR), rS),       rS)
    end
    if return_nFs   
        return nFs_Au, nFs_YSZ, nFs_tot
    else
        return poisson_eq
    end
end


























































function HALF_test_analytical(data, V_tot, V_ISR; verbose=true)

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
    #g = HALF_test_analytical(data, V_tot, 666, verbose=false)
    g = HALF_test_analytical_AYA_Sl(data, V_tot, 666, verbose=false)
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


function HALF_get_analytic_V_ISR__AYA_Sl(data, V_tot)
    g = HALF_test_analytical_AYA_Sl(data, V_tot, 666, verbose=false)
    to_optimize(x) = g(x[1])^2
    Optim.optimize(to_optimize, -10, 10, GoldenSection()).minimizer
end

struct TestError <: Exception
    text::String
end

Base.showerror(io::IO, e::TestError) = print(io, e.text)

function test_one_AYA_Sl_cell(STcell::AYA_Sl)
    data = STcell.system.physics.data
    testing_id = STcell.system.grid.components[ExtendableGrids.BFaceNodes][3]
    phi_ISR = STcell.U[1, testing_id]
    phi_L = STcell.U[1, 1]
    phi_R = STcell.U[1, Int(round(end/2 - 2))]
    V_tot = phi_L - phi_R
    V_ISR = phi_ISR - phi_R
    #@show phi_R
    rel_phi_err = (HALF_get_analytic_V_ISR__AYA_Sl(data, V_tot) - V_ISR)/V_ISR
    return rel_phi_err
end

function test_AYA_Sl_cell()
    STcell = AYA_Sl()
    stationary_update!(STcell,Dict([
				:AYSZ => 0.0,
                :AYSZs => 0.0,
				
				:GA => 0.2*e0,
				:Ge => 0.0*e0,
				:alpha => 0.01,
                :alphas => 0.01,

                :DYSZ => 1.0e-8,
				:kA => 1.0e25,
				
				#:kA => 0.0,
				#:boundary_charge_fac => 0.0
			]))

    function perform_it(bias, STcell)
        if !test_partial_charges(STcell)
            throw(TestError("TestError: Partial charges are not allright"))
        end
        rel_phi_err = test_one_AYA_Sl_cell(STcell::AYA_Sl)
        #@show bias, rel_phi_err
        push!(rel_err_storage, rel_phi_err)
    end

    rel_err_storage = []
    biasshow!(STcell, bend=0.2, bstep=0.01, tend=1e-2, callback=perform_it)


    return rel_err_storage
end