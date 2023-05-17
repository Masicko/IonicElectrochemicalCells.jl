
#=

Asymetric variant for density of bulk quantities

=#
@composite mutable struct AYA_Sl <: AbstractCell
    CommonCell...
end

#function AYA_Sl(grid=cell1D(electrode_thickness, electrolyte_thickness, electrode_thickness, dmin=1e-16))
function AYA_Sl(grid=doublehalf_cell1D(electrode_thickness, electrolyte_thickness, electrode_thickness, dmin=1e-16))
    new = AYA_Sl(0.0, 0.0, 0.0, 0.0)
    new.parameters = Model((SharedParams(), YSZparams(), Auparams(), ISRparameters()))
    new.system = AYA_Sl_system(grid, AueDensity=e_BoltzmannAu_ne)
    new.U = inival(new)
    new.Uold = copy(new.U)
    return new
end

e_BoltzmannAu_ne(data, V, S) = enLAu(S) * exp(-ze * e0 / kB / data.T * V)

AYA_Sl_physics = function (; AueDensity=e_BoltzmannAu_ne)
    return VoronoiFVM.Physics(
        data=AYALG1iBoltzmannData(),
        
        flux=function (f, u, edge, data)
            if edge.region in (Ω_Aul, Ω_Aur)
                f[ipsi] = Fick_flux(u[ipsi, 1], u[ipsi, 2], epsAu)
            elseif edge.region in (e_Ω_YSZl, e_Ω_YSZr)
                f[iyV] = LGS_flux(u[iyV, 1], u[iyV, 2], u[ipsi, 1], u[ipsi, 2], data.DYSZ*get_S(edge, data), e0 * zV / kB / data.T, data.AYSZ * e0 / kB / data.T) # interaction energy A enters here
                f[ipsi] = Fick_flux(u[ipsi, 1], u[ipsi, 2], epsYSZ)
            end
        end,
        #
        reaction=function (f, u, node, data)
            if node.region == Ω_Aul
                f[ipsi] = -e_Au_charge_density(AueDensity(data, u[ipsi]- data.bias, data.SL), data.SL)
            elseif node.region in (e_Ω_YSZl, e_Ω_YSZr)
                f[ipsi] = -e_YSZ_charge_density(e_nVmax(data.alpha, get_S(node, data)) * u[iyV], get_S(node, data)) 
            elseif node.region == Ω_Aur
                f[ipsi] = -e_Au_charge_density(AueDensity(data, u[ipsi] - 0, data.SR), data.SR)
            end
        end,
        #
        storage=function (f, u, node, data)
            f[ipsi] = 0.0
            if node.region in (e_Ω_YSZl, e_Ω_YSZr)
                f[iyV] = u[iyV]*get_S(node, data)
            end
        end,

        breaction=function (f, u, bnode, data)
            #if bnode.region in e_ISR
            if bnode.region in ISR
                # !! notation: ν is an outer normal vector
                # Adsorption of vacancies 
                KVsq = exp(-data.GA / data.T / kB / 2 + (data.AYSZ * u[iyV] - data.AYSZs * u[iyVs]) / 2) # sqrt(KV)
                ReducedRateA = KVsq * u[iyV] * (1.0 - u[iyVs]) - 1 / KVsq * u[iyVs] * (1.0 - u[iyV])
                # INFO the reaction rate is calculated in [# vacancies/cross section area]
                # and its contribution to both coverages isn't thus far scaled appropriately
                # However, this is no problem for a blocking electrode in equilibrium...
                # the fix should look something like
                # ... a check is needed though
                Sl = get_S(bnode, data)
                # implementation for bspecies
                # bstorage + breaction = 0
                f[iyVs] = -data.kA * ReducedRateA / nVmaxs(data.alphas, Sl)
                # implementation for species
                # - j ̇ν + breaction = 0
                f[iyV] = Sl*data.kA * ReducedRateA / e_nVmax(data.alpha, Sl)
                # equilibrium of electrons
                # V = (bnode.region == Γ_YSZl ? data.bias : 0.0) - u[ipsi]
                # V = reduced_voltage(u, bnode, data)
                # TODO add the difference of the electrostatic potential or its derivative*thickness of the ISR to the "equilibrium constant for electrons"
                # FIXME implicitly assuming the Boltzmann statistics for the ISR electrons
                # nes = nLAus(S) / nLAu * exp(-data.Ge / kB / data.T) * AueDensity(data, V) # [# electrons/ ISR area]
                nes = e_ISR_electrondensity(u, bnode, data) # [# electrons/ ISR area]
                # TODO ISRthickness = (bnode.region  == Γ_YSZl ? data.dL : data.dR)
                ###
                # The surface Poisson equation is consistent with [BSE2018]
                # (-ε_+ ∇ψ_+)⋅ν_+ (-ε_- ∇ψ_-)⋅ν_- + ISR_chargedensity = 0
                # However, the equation for normal fluxes and the breaction for an internal node between two REVs is implemented as
                # - (j_+ ̇ν_+ + j_- ̇ν_-) + breaction = 0
                # thus ISR_chargedensity below, correctly, enters with a negative sign !!!
                f[ipsi] = (data.boundary_charge_fac)*-ISR_chargedensity(nes, nVmaxs(data.alphas, Sl ) * u[iyVs], Sl)
                #
            end
        end,
        #
        bstorage=function (f, u, bnode, data)
            if bnode.region in ISR
                    f[iyVs] = u[iyVs]
            end
        end
    )
end

# AYA_Sl_system constructor
function AYA_Sl_system(grid; AueDensity=e_BoltzmannAu_ne)
    system = VoronoiFVM.System(grid, AYA_Sl_physics(AueDensity=AueDensity), unknown_storage=:sparse)
    #
    enable_species!(system, ipsi, [Ω_Aul, e_Ω_YSZl, e_Ω_YSZr, Ω_Aur])
    enable_species!(system, iyV, [e_Ω_YSZl, e_Ω_YSZr])
    enable_boundary_species!(system, iyVs, [Γ_YSZl, Γ_YSZr])
    #
    boundary_dirichlet!(system, ipsi, Γ_Aul, 0.0)
    boundary_dirichlet!(system, ipsi, Γ_Aur, 0.0)
    #
    #check_allocs!(system, true)
    return system
end

# initial values for AYA_Sl_system
function AYA_Sl_inival(system)
    inival = VoronoiFVM.unknowns(system, inival=e_yV_en(system.physics.data.alpha, 1.0))
    inival[ipsi, :] .= 0.0
    boundary_ids = system.grid.components[ExtendableGrids.BFaceNodes]
    ISR_L_id = boundary_ids[3]
	inival[iyVs, ISR_L_id] = yVs_en(system.physics.data.alphas, system.physics.data.SL)
    if length(boundary_ids) > 3
        ISR_R_id = boundary_ids[4]
        inival[iyVs, ISR_R_id] = yVs_en(system.physics.data.alphas, system.physics.data.SR)
    end
    return inival
end

inival(cell::AYA_Sl) = AYA_Sl_inival(cell.system)









#################################################################################
######################### postprocessing stuff ##################################


mutable struct dummy_bnode
	region
end

function get_partial_charges(cell::AYA_Sl; only_DL=false)
    data = cell.system.physics.data
    BFNodes = cell.system.grid.components[ExtendableGrids.BFaceNodes]
    b1_id = BFNodes[3]
    b2_id = BFNodes[4]

    U = cell.U

    partial_charge_dict = Dict()
    
    partial_charge_dict[:bQ_YSZ_L] = YSZ_surface_charge(U[3, b1_id], data.alphas, data.SL)
    partial_charge_dict[:bQ_Au_L] = Au_surface_charge(
		(1/nLAus(data.SL))*ISR_electrondensity(
			U[:, b1_id],  dummy_bnode(Γ_YSZl), data 
		), 
		data.SL
	)

    partial_charge_dict[:bQ_YSZ_R] = YSZ_surface_charge(U[3, b2_id], data.alphas, data.SR)
    partial_charge_dict[:bQ_Au_R] = Au_surface_charge(
		(1/nLAus(data.SR))*ISR_electrondensity(
			U[:, b2_id],  dummy_bnode(Γ_YSZr), data 
		), 
		data.SR
	)

    bulk_symbols =      [:nF_Au_L, :nF_YSZ_L,   :nF_YSZ_R,    :nF_Au_R]
    bulk_identifiers =  [Ω_Aul,     e_Ω_YSZl,   e_Ω_YSZr,     Ω_Aur]
    bulk_charges = bulk_charge(cell)[1, :]

    [partial_charge_dict[bulk_symbols[i]] = bulk_charges[bulk_identifiers[i]] for i in 1:4]  
    
    if only_DL
        return (
            partial_charge_dict[:nF_Au_L] + partial_charge_dict[:bQ_Au_L],
            partial_charge_dict[:nF_YSZ_L] + partial_charge_dict[:bQ_YSZ_L],
            partial_charge_dict[:nF_YSZ_R] + partial_charge_dict[:bQ_YSZ_R],
            partial_charge_dict[:nF_Au_R] + partial_charge_dict[:bQ_Au_R]
        )
    end

    return partial_charge_dict
end

function test_partial_charges_full_cell(pch)
    
    tol = 1.0e-11
    return (
        # DLs
        abs(pch[:nF_Au_L] + pch[:bQ_Au_L] + pch[:nF_YSZ_L] + pch[:bQ_YSZ_L]) < tol 
        && abs(pch[:nF_Au_R] + pch[:bQ_Au_R] + pch[:nF_YSZ_R] + pch[:bQ_YSZ_R]) < tol
        
        # L ans R incoming charge from electrons
        && abs(pch[:nF_Au_L] + pch[:bQ_Au_L] + pch[:nF_Au_R] + pch[:bQ_Au_R]) < tol

        # YSZ
        && abs(pch[:nF_YSZ_L] + pch[:bQ_YSZ_L] + pch[:nF_YSZ_R] + pch[:bQ_YSZ_R]) < tol
    )
end

function test_partial_charges(cell::AYA_Sl)
    pch = get_partial_charges(cell)
    test_partial_charges_full_cell(pch)
end


function get_stored_charge(cell::AYA_Sl)
    # 2* because get_charge returns only Au side of left DL and there is also equivalent charge on Au on right DL.
    2*get_half_DL_charge(cell, "l")
end

function get_half_DL_charge(cell::AYA_Sl, side="l")
    QrAuL = VoronoiFVM.integrate(cell.system, cell.system.physics.reaction, cell.U) # nspec x nregion
    
    if side == "l"
        bulk, bnd = Ω_Aul, Γ_YSZl  # Γ_YSZl
    elseif side == "r"
        bulk, bnd = Ω_Aur, Γ_YSZr  # Γ_YSZr
    else
        println("side = $(side) which is not \"r\" or \"l\"")
        return throw(Exception())
    end
    
    bnd_index = cell.system.grid.components[BFaceNodes][bnd]
    data = cell.system.physics.data
    
    S = (bnd == Γ_YSZl ? data.SL : data.SR)
    V = cell.U[ipsi, bnd_index] - (bnd == Γ_YSZl ? cell.U[ipsi, 1] : 0.0)

    nes = nLAus(S) / enLAu(S) * exp(-data.Ge / kB / data.T) * e_BoltzmannAu_ne(data, V, S) # [# electrons/ ISR area]
    
    ISRcontribution = (data.boundary_charge_fac)*e0*(nLAus(S) - nes)
    #@show V, ISRcontribution, -QrAuL[ipsi, bulk]

    return -QrAuL[ipsi, bulk] + ISRcontribution # the sign of the surface charge is opposite, see surface Poisson
end