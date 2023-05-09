
# fluxes used by all YSZ models
function flux_test(sys)
    edge = nothing
    data = sys.physics.data
    flux_wrapper!(f, u) = sys.physics.flux!(f, u, edge, data)
    u_kl = [1.0 2.0; 3.0 2.0; 2.0 1.0]
    f_k = [0.8, 0.8, 0.8]
    @show ForwardDiff.jacobian(flux_wrapper!, f_k, u_kl)
end        

function YSZfluxes!(f, u, edge, data)
    f[iyV] = LGS_flux(u[iyV, 1], u[iyV, 2], u[ipsi, 1], u[ipsi, 2], data.DYSZ, e0 * zV / kB / data.T, data.AYSZ * e0 / kB / data.T) # interaction energy A enters here
    f[ipsi] = Fick_flux(u[ipsi, 1], u[ipsi, 2], epsYSZ)
end

function YSZstorage!(f, u, node, data)
    f[iyV] = u[iyV]
    f[ipsi] = 0.0
end

function YSZreaction!(f, u, node, data)
    f[iyV] = 0.0
    f[ipsi] = -YSZ_charge_density(nVmax(data.alpha) * u[iyV])
end





































# @composite mutable struct YHElectrolyte <: AbstractCell
#     CommonCell...
# end
# function YHElectrolyte(grid=half_electrolyte_1D(electrolyte_length, dmin))
#     new = YHElectrolyte(Model((SharedParams, YSZparams)), bulk_lattice_half_cell(grid=grid), 0.0, 0.0)
#     new.U = inival(new)
#     new.Uold = copy(new.U)
#     return new
# end

# parameters2VoronoiData(Model((SharedParams(), YSZparams())), :parametersY)()

# function bulk_lattice_half_cell(; grid=nothing)
#     grid = half_electrolyte_1D(electrolyte_length, dmin)  # TODO verify the indexing of the grid
#     data = parametersY()
#     physics = VoronoiFVM.Physics(
#         data=data,
#         num_species=2,
#         flux=YSZfluxes!,
#         #
#         storage=YSZstorage!,
#         ##
#         reaction=YSZreaction!
#     )
#     sys = VoronoiFVM.System(grid, physics, unknown_storage=:sparse)


#     enable_species!(sys, iyV, [Ω_YSZ])
#     enable_species!(sys, ipsi, [Ω_YSZ])
#     # electroneutral Dirichlet and zero Neumann for iyV
#     boundary_dirichlet!(sys, iyV, 2, electroneutral_nV / nVmax(data.alpha))
#     # Dirichlet BCs for ipsi
#     boundary_dirichlet!(sys, ipsi, 1, data.bias)
#     boundary_dirichlet!(sys, ipsi, 2, 0.0)

#     return sys
# end

# function inival(cell::YHElectrolyte)
#     inival = VoronoiFVM.unknowns(cell.system)
#     inival[ipsi, :] .= 0.0
#     inival[iyV, :] .= electroneutral_nV / nLYSZ
#     return inival
# end


# function current_integrator(sys, U)
#     factory = VoronoiFVM.TestFunctionFactory(sys)
#     tfL = testfunction(factory, 2, 1)
#     return VoronoiFVM.integrate_stdy(sys, tfL, U)
# end

# function impedance_current(sys, U)
#     return current_integrator(sys, U)[ipsi]
# end



































































#=

The full transient gPNP problem for the blocking Au|YSZ|Au cell

=#
@composite mutable struct AYABoltzmann <: AbstractCell
    CommonCell...
end

parameters2VoronoiData(Model((SharedParams(), YSZparams(), Auparams())), :AYABoltzmannData)()

function AYABoltzmann(grid=cell1D(electrode_thickness, electrolyte_thickness, electrode_thickness, dmin=1e-12))
    new = AYABoltzmann(Model((SharedParams(), YSZparams(), Auparams())), ayasystem(grid)[1], 0.0, 0.0)
    new.U = inival(new)
    new.Uold = copy(new.U)
    return new
end


function ayasystem(grid)
    system = VoronoiFVM.System(grid, unknown_storage=:sparse)
    psi = VoronoiFVM.ContinuousQuantity(system, collect(bulk_domains))
    dspec = VoronoiFVM.DiscontinuousQuantity(system, collect(bulk_domains))
    function flux(f, u, edge, data)
        if edge.region in (Ω_Aul, Ω_Aur)
            f[psi] = Fick_flux(u[psi, 1], u[psi, 2], epsAu)
            f[dspec] = PB_flux(u[dspec, 1], u[dspec, 2], u[psi, 1], u[psi, 2], data.DAu, e0 * ze / kB / data.T)
        elseif edge.region == Ω_YSZ
            #f[psi] = Fick_flux(u[psi, 1], u[psi, 2], epsYSZ)
            #f[dspec] = LGS_flux(u[dspec, 1], u[dspec, 2], u[psi, 1], u[psi, 2], data.DYSZ, e0 * zV / kB / data.T, data.AYSZ * e0 / kB / data.T)
            YSZfluxes!(f,u,edge,data)
        end
    end

    function storage(f, u, node, data)
        f[psi] = 0.0
        f[dspec] = u[dspec]
    end

    function reaction(f, u, node, data)
        if node.region in (Ω_Aul, Ω_Aur)
            f[psi] = -Au_charge_density(nLAu * u[dspec])
        elseif node.region == Ω_YSZ
            f[psi] = -YSZ_charge_density(nVmax(data.alpha) * u[dspec])
        end
    end

    VoronoiFVM.physics!(system,
        VoronoiFVM.Physics(
            # data=parameters(),
            data=AYABoltzmannData(),
            flux=flux,
            reaction=reaction,
            storage=storage,
        )
    )
    boundary_dirichlet!(system, psi, 1, 0.0)
    boundary_dirichlet!(system, psi, 2, 0.0)
    #
    boundary_dirichlet!(system, dspec, 1, neutral_yeAU)
    boundary_dirichlet!(system, dspec, 2, neutral_yeAU)
    #
    check_allocs!(system, true)
    return system, psi, dspec
end

#AYAsystem() = ayasystem()[1]

function aya_inival(system)
    inival = VoronoiFVM.unknowns(system, inival=nVmax(0.0) / nVmax(system))
    inival[1, :] .= 0.0
    if size(inival)[1] > 2
        inival[2, :] .= 1.0
        inival[4, :] .= 1.0
    end
    return inival
end
inival(cell::AYABoltzmann) = aya_inival(cell.system)







































#=
Reduced problem: 
- YSZ -- transient gPNP for ionic vacancies
- Au -- nonlinear Poisson problem

=#
# FIXME solve the math in #= ... =# comments
@composite mutable struct AYALGBoltzmann <: AbstractCell
    CommonCell...
end

function AYALGBoltzmann(grid=cell1D(electrode_thickness, electrolyte_thickness, electrode_thickness, dmin=1e-16))
    new = AYALGBoltzmann(Model((SharedParams(), YSZparams(), Auparams())), ayasystemLG(grid, AueDensity=BoltzmannAu_ne), 0.0, 0.0)
    new.U = inival(new)
    new.Uold = copy(new.U)
    return new
end


ayaLGphys = function (; AueDensity=BoltzmannAu_ne)
    return VoronoiFVM.Physics(
        # data=parameters(),
        data=AYABoltzmannData(),
        #
        flux=function (f, u, edge, data)
            if edge.region in (Ω_Aul, Ω_Aur)
                f[ipsi] = Fick_flux(u[ipsi, 1], u[ipsi, 2], epsAu)
            elseif edge.region == Ω_YSZ
                YSZfluxes!(f,u,edge,data)
            end
        end,
        #   
        reaction=function (f, u, node, data)
            if node.region == Ω_Aul
                f[ipsi] = -Au_charge_density(AueDensity(data, u[ipsi]- data.bias))
            elseif node.region == Ω_YSZ
                f[ipsi] = -YSZ_charge_density(nVmax(data.alpha) * u[iyV])   
            elseif node.region == Ω_Aur
                f[ipsi] = -Au_charge_density(AueDensity(data, u[ipsi] - 0))
            end
        end,
        #
        storage=function (f, u, node, data)
            f[ipsi] = 0.0
            if node.region == 2
                #f[iyV] = u[iyV]
                YSZstorage!(f, u, node, data)
            end
        end,
    )
end

# V = phi - phi_B

# Thomas-Fermi-Dirac model of electrons in metal 
# TODO check -> correct -> test
# const Ck = 3*Planck_constant^2/40/mₑ/a0^2*(3/pi)^(2/3) 
# TFDAue(data, V) = (nLAu^(2/3) + V/Ck)^(3/2) 

function ayasystemLG(grid; AueDensity=BoltzmannAu_ne)
    system = VoronoiFVM.System(grid, ayaLGphys(AueDensity=AueDensity), unknown_storage=:sparse)
    #
    enable_species!(system, ipsi, collect(bulk_domains))
    enable_species!(system, iyV, [Ω_YSZ])
    #
    boundary_dirichlet!(system, ipsi, 1, 0.0)
    boundary_dirichlet!(system, ipsi, 2, 0.0)
    #
    check_allocs!(system, true)
    return system
end

function ayaLGinival(system)
    inival = VoronoiFVM.unknowns(system, inival=nVmax(0.0) / nVmax(system))
    inival[ipsi, :] .= 0.0
    return inival
end

inival(cell::AYALGBoltzmann) = ayaLGinival(cell.system)

function get_charge(cell::AYALGBoltzmann, side)
    QrAuL = VoronoiFVM.integrate(cell.system, cell.system.physics.reaction, cell.U) # nspec x nregion
    
    if side == "l"
        bulk, bnd = Ω_Aul, Γ_YSZl  # Γ_YSZl
    elseif side == "r"
        bulk, bnd = Ω_Aur, Γ_YSZr  # Γ_YSZr
    else
        println("side = $(side) which is not \"r\" or \"l\"")
        return throw(Exception())
    end
    
    #bnd_index = cell.system.grid.components[BFaceNodes][bnd]
    #data = cell.system.physics.data
    
    #V = (bnd == Γ_YSZl ? cell.U[ipsi, 1] : 0.0) - cell.U[ipsi, bnd_index]
    #nes = nLAus(S) / nLAu * exp(-data.Ge / kB / data.T) * BoltzmannAu_ne(data, V) # [# electrons/ ISR area]
    #ISRcontribution = e0*nes

    #@show V, ISRcontribution, -QrAuL[ipsi, bulk]

    return -QrAuL[ipsi, bulk] # the sign of the surface charge is opposite, see surface Poisson
end














































#=

the previous case extended with interface species
chemical potential of the surface vacancies

=#
# TODO governing eqns for the interface species
@composite mutable struct AYALG1iBoltzmann <: AbstractCell
    CommonCell...
end


function AYALG1iBoltzmann(grid=cell1D(electrode_thickness, electrolyte_thickness, electrode_thickness, dmin=1e-16))
    new = AYALG1iBoltzmann(0.0, 0.0, 0.0, 0.0)
    new.parameters = Model((SharedParams(), YSZparams(), Auparams(), ISRparameters()))
    new.system = ayasystemLG1i(grid, AueDensity=BoltzmannAu_ne)
    new.U = inival(new)
    new.Uold = copy(new.U)
    return new
end

parameters2VoronoiData(Model((SharedParams(), YSZparams(), Auparams(), ISRparameters())), :AYALG1iBoltzmannData)()

ayaLGphys1i = function (; AueDensity=BoltzmannAu_ne)
    ayaLGp = ayaLGphys(AueDensity=AueDensity)
    return VoronoiFVM.Physics(
        data=AYALG1iBoltzmannData(),
        # bulk properties inherited from ayaLGphysics
        flux=ayaLGp.flux,
        reaction=ayaLGp.reaction,
        storage=ayaLGp.storage,
        breaction=function (f, u, bnode, data)
            if bnode.region in ISR
                # !! notation: ν is an outer normal vector
                # Adsorption of vacancies 
                # TODO add the difference of the electrostatic potential or its derivative*thickness of the ISR to the "equilibrium constant for vacancies"
                KVsq = exp(-data.GA / data.T / kB / 2 + (data.AYSZ * u[iyV] - data.AYSZs * u[iyVs]) / 2) # sqrt(KV)
                ReducedRateA = KVsq * u[iyV] * (1.0 - u[iyVs]) - 1 / KVsq * u[iyVs] * (1.0 - u[iyV])
                # INFO the reaction rate is calculated in [# vacancies/cross section area]
                # and its contribution to both coverages isn't thus far scaled appropriately
                # However, this is no problem for a blocking electrode in equilibrium...
                # the fix should look something like
                # ... a check is needed though
                # AreaRatio = (bnode.region == Γ_YSZl ? data.SL : data.SR)
                AreaRatio = ISR_arearatio(bnode, data)
                # implementation for bspecies 
                # bstorage + breaction = 0
                f[iyVs] = -data.kA * ReducedRateA / nVmaxs(data.alphas, AreaRatio)
                # implementation for species
                # - j ̇ν + breaction = 0
                f[iyV] = data.kA * ReducedRateA / nVmax(data.alpha)
                # equilibrium of electrons
                # V = (bnode.region == Γ_YSZl ? data.bias : 0.0) - u[ipsi]
                # V = reduced_voltage(u, bnode, data)
                # TODO add the difference of the electrostatic potential or its derivative*thickness of the ISR to the "equilibrium constant for electrons"
                # FIXME implicitly assuming the Boltzmann statistics for the ISR electrons
                # nes = nLAus(S) / nLAu * exp(-data.Ge / kB / data.T) * AueDensity(data, V) # [# electrons/ ISR area]
                nes = ISR_electrondensity(u, bnode, data) # [# electrons/ ISR area]
                # TODO ISRthickness = (bnode.region  == Γ_YSZl ? data.dL : data.dR)
                ###
                # The surface Poisson equation is consistent with [BSE2018]
                # (-ε_+ ∇ψ_+)⋅ν_+ (-ε_- ∇ψ_-)⋅ν_- + ISR_chargedensity = 0
                # However, the equation for normal fluxes and the breaction for an internal node between two REVs is implemented as
                # - (j_+ ̇ν_+ + j_- ̇ν_-) + breaction = 0
                # thus ISR_chargedensity below, correctly, enters with a negative sign !!!
                f[ipsi] = (data.boundary_charge_fac)*-ISR_chargedensity(nes, nVmaxs(data.alphas, AreaRatio ) * u[iyVs], AreaRatio)
                #
            end
        end,
        #
        bstorage=function (f, u, bnode, data)
            if bnode.region in ISR #[Γ_YSZl,Γ_YSZr]
                f[iyVs] = u[iyVs]
            end
        end
    )
end

# ayasystemLGi1 constructor
function ayasystemLG1i(grid; AueDensity=BoltzmannAu_ne)
    system = VoronoiFVM.System(grid, ayaLGphys1i(AueDensity=AueDensity), unknown_storage=:sparse)
    #
    enable_species!(system, ipsi, collect(bulk_domains))
    enable_species!(system, iyV, [Ω_YSZ])
    enable_boundary_species!(system, iyVs, [Γ_YSZl, Γ_YSZr])
    #
    boundary_dirichlet!(system, ipsi, 1, 0.0)
    boundary_dirichlet!(system, ipsi, 2, 0.0)
    #
    #check_allocs!(system, true)
    return system
end
# initial values for ayasystemLGi1
function ayaLG1i_inival(system)
    inival = VoronoiFVM.unknowns(system, inival=nVmax(0.0) / nVmax(system))
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

inival(cell::AYALG1iBoltzmann) = ayaLG1i_inival(cell.system)

mutable struct whatever
    region
end

function get_partial_charges(cell::AYALG1iBoltzmann)
    data = cell.system.physics.data
    X = cell.system.grid.components[ExtendableGrids.Coordinates]
    BFNodes = cell.system.grid.components[ExtendableGrids.BFaceNodes]
    b1_id = BFNodes[3]
    b2_id = BFNodes[4]

    U = cell.U

    partial_charge_dict = Dict()

    # surface
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

    # bulk
    get_charge_from_sys(yV) = YSZ_charge_density(
		nVmax(data.alpha) * yV
	)
    charge_list = get_charge_from_sys.(U[2, b1_id : b2_id])
    X_list = X[b1_id : b2_id]
    center = Int(round(length(X_list)/2))
    partial_charge_dict[:nF_YSZ_L] = integrate_arrays_by_trapezoid(X_list[1:center], charge_list[1:center])
    partial_charge_dict[:nF_YSZ_R] = integrate_arrays_by_trapezoid(X_list[center:end], charge_list[center:end])

    bulk_symbols =      [:nF_Au_L,  :nF_Au_R]
    bulk_identifiers =  [Ω_Aul,     Ω_Aur]
    bulk_charges = bulk_charge(cell)[1, :]

    [partial_charge_dict[bulk_symbols[i]] = bulk_charges[bulk_identifiers[i]] for i in 1:length(bulk_identifiers)]  
    
    return partial_charge_dict
end


















#=

HALF cell with everything... but also with analytical solution for the equilibrium case
=#
@composite mutable struct AYALG1iBoltzmann_HALF <: AbstractCell
    CommonCell...
end

function AYALG1iBoltzmann_HALF(grid=half_cell1D(electrode_thickness, electrolyte_thickness, electrode_thickness, dmin=1.0e-16))
    new = AYALG1iBoltzmann_HALF(0.0, 0.0, 0.0, 0.0)
    new.parameters = Model((SharedParams(), YSZparams(), Auparams(), ISRparameters()))
    new.system = ayasystemLG1i_HALF(grid)
    new.U = inival(new)
    new.Uold = copy(new.U)
    return new
end

function ayasystemLG1i_HALF(grid, AueDensity=BoltzmannAu_ne)
    system = VoronoiFVM.System(grid, ayaLGphys1i(AueDensity=AueDensity), unknown_storage=:sparse)
    #
    enable_species!(system, ipsi, [Ω_Aul, Ω_YSZ])
    enable_species!(system, iyV, [Ω_YSZ])
    enable_boundary_species!(system, iyVs, [Γ_YSZl])
    #
    boundary_dirichlet!(system, ipsi, 1, 0.0)
    boundary_dirichlet!(system, ipsi, 2, 0.0)
    #boundary_dirichlet!(system, iyV, 2, yV_en(system.physics.data.alpha))
    #
    check_allocs!(system, true)
    return system
end

inival(cell::AYALG1iBoltzmann_HALF) = ayaLG1i_inival(cell.system)



















































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

function get_S(edge, data)
    if edge.region in [Ω_Aul, e_Ω_YSZl]
        return data.SL
    elseif edge.region in [Ω_Aur, e_Ω_YSZr]
        return data.SR
    end
end

e_BoltzmannAu_ne(data, V, S) = enLAu(S) * exp(-ze * e0 / kB / data.T * V)

AYA_Sl_physics = function (; AueDensity=e_BoltzmannAu_ne)
    return VoronoiFVM.Physics(
        data=AYALG1iBoltzmannData(),
        # bulk properties inherited from ayaLGphysics
        #flux=ayaLGp.flux,
        #reaaction=ayaLGp.reaction,
        #storage=ayaLGp.storage,
        
        flux=function (f, u, edge, data)
            if edge.region in [Ω_Aul, Ω_Aur]
                f[ipsi] = Fick_flux(u[ipsi, 1], u[ipsi, 2], epsAu)
            elseif edge.region in [e_Ω_YSZl, e_Ω_YSZr]
                f[iyV] = LGS_flux(u[iyV, 1], u[iyV, 2], u[ipsi, 1], u[ipsi, 2], data.DYSZ*get_S(edge, data), e0 * zV / kB / data.T, data.AYSZ * e0 / kB / data.T) # interaction energy A enters here
                f[ipsi] = Fick_flux(u[ipsi, 1], u[ipsi, 2], epsYSZ)
            end
        end,
        #   
        reaction=function (f, u, node, data)
            if node.region == Ω_Aul
                #f[ipsi] = -Au_charge_density(AueDensity(data, u[ipsi]- data.bias))
                f[ipsi] = -e_Au_charge_density(AueDensity(data, u[ipsi]- data.bias, data.SL), data.SL)
            elseif node.region in [e_Ω_YSZl, e_Ω_YSZr] #[eΩ_YSZl, eΩ_YSZr]
                #f[ipsi] = -YSZ_charge_density(nVmax(data.alpha) * u[iyV])  
                f[ipsi] = -e_YSZ_charge_density(e_nVmax(data.alpha, get_S(node, data)) * u[iyV], get_S(node, data)) 
            elseif node.region == Ω_Aur
                #f[ipsi] = -Au_charge_density(AueDensity(data, u[ipsi] - 0))
                f[ipsi] = -e_Au_charge_density(AueDensity(data, u[ipsi] - 0, data.SR), data.SR)
            end
        end,
        #
        storage=function (f, u, node, data)
            f[ipsi] = 0.0
            if node.region in [e_Ω_YSZl, e_Ω_YSZr]
                f[iyV] = u[iyV]*get_S(node, data)
            end
        end,

        # flux=function (f, u, edge, data)
        #     if edge.region in (eΩ_Aul, eΩ_Aur)
        #         f[ipsi] = Fick_flux(u[ipsi, 1], u[ipsi, 2], epsAu)
        #     elseif edge.region in (eΩ_YSZl, eΩ_YSZr)
        #         f[iyV] = LGS_flux(u[iyV, 1], u[iyV, 2], u[ipsi, 1], u[ipsi, 2], data.DYSZ*get_S(edge, data), e0 * zV / kB / data.T, data.AYSZ * e0 / kB / data.T) # interaction energy A enters here
        #         f[ipsi] = Fick_flux(u[ipsi, 1], u[ipsi, 2], epsYSZ)
        #     end
        # end,
        # reaction=function (f, u, node, data)
        #     if node.region == eΩ_Aul
        #         f[ipsi] = -e_Au_charge_density(AueDensity(data, u[ipsi]- data.bias, data.SL), data.SL)
        #     elseif node.region == eΩ_YSZl
        #         f[ipsi] = -e_YSZ_charge_density(e_nVmax(data.alpha, data.SL) * u[iyV], data.SL)
        #     elseif node.region == eΩ_YSZr
        #         f[ipsi] = -e_YSZ_charge_density(e_nVmax(data.alpha, data.SR) * u[iyV], data.SR)   
        #     elseif node.region == eΩ_Aur
        #         f[ipsi] = -e_Au_charge_density(AueDensity(data, u[ipsi] - 0, data.SR), data.SR)
        #     end
        # end,
        # storage=function (f, u, node, data)
        #     f[ipsi] = 0.0
        #     if node.region in (eΩ_YSZl, eΩ_YSZr)
        #         f[iyV] = u[iyV]*get_S(node, data)
        #     end
        # end,
        breaction=function (f, u, bnode, data)
            #if bnode.region in e_ISR
            if bnode.region in ISR
                    # !! notation: ν is an outer normal vector
                # Adsorption of vacancies 
                # TODO add the difference of the electrostatic potential or its derivative*thickness of the ISR to the "equilibrium constant for vacancies"
                KVsq = exp(-data.GA / data.T / kB / 2 + (data.AYSZ * u[iyV] - data.AYSZs * u[iyVs]) / 2) # sqrt(KV)
                ReducedRateA = KVsq * u[iyV] * (1.0 - u[iyVs]) - 1 / KVsq * u[iyVs] * (1.0 - u[iyV])
                # INFO the reaction rate is calculated in [# vacancies/cross section area]
                # and its contribution to both coverages isn't thus far scaled appropriately
                # However, this is no problem for a blocking electrode in equilibrium...
                # the fix should look something like
                # ... a check is needed though
                # AreaRatio = (bnode.region == Γ_YSZl ? data.SL : data.SR)
                AreaRatio = ISR_arearatio(bnode, data)
                # implementation for bspecies
                # bstorage + breaction = 0
                f[iyVs] = -data.kA * ReducedRateA / nVmaxs(data.alphas, AreaRatio)
                # implementation for species
                # - j ̇ν + breaction = 0
                f[iyV] = data.kA * ReducedRateA / e_nVmax(data.alpha, AreaRatio)
            #f[iyV] = data.kA * ReducedRateA / nVmax(data.alpha)
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
                f[ipsi] = (data.boundary_charge_fac)*-ISR_chargedensity(nes, nVmaxs(data.alphas, AreaRatio ) * u[iyVs], AreaRatio)
                #
            end
        end,
        #
        bstorage=function (f, u, bnode, data)
            #if bnode.region in e_ISR #[Γ_YSZl,Γ_YSZr]
            if bnode.region in ISR #[Γ_YSZl,Γ_YSZr]
                    f[iyVs] = u[iyVs]
            end
        end
    )
end

# AYA_Sl_system constructor
function AYA_Sl_system(grid; AueDensity=e_BoltzmannAu_ne)
    system = VoronoiFVM.System(grid, AYA_Sl_physics(AueDensity=AueDensity), unknown_storage=:sparse)
    #
    #enable_species!(system, ipsi, collect([eΩ_Aul, 5, eΩ_Aur]))
    #enable_species!(system, iyV, [5])
    #enable_boundary_species!(system, iyVs, [eΓ_YSZl, eΓ_YSZr])
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

mutable struct dummy_bnode
	region
end

function get_partial_charges(cell::AYA_Sl)
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
    
    return partial_charge_dict
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

    nes = nLAus(S) / nLAu(S) * exp(-data.Ge / kB / data.T) * e_BoltzmannAu_ne(data, V, S) # [# electrons/ ISR area]
    
    ISRcontribution = (data.boundary_charge_fac)*e0*(nLAus(S) - nes)
    #@show V, ISRcontribution, -QrAuL[ipsi, bulk]

    return -QrAuL[ipsi, bulk] + ISRcontribution # the sign of the surface charge is opposite, see surface Poisson
end



















function get_stored_charge(cell::AYALG1iBoltzmann)
    # 2* because get_charge returns only Au side of left DL and there is also equivalent charge on Au on right DL.
    2*get_half_DL_charge(cell, "l")
end

function AuL_charge(cell::AYALG1iBoltzmann)
    get_half_DL_charge(cell, "l")
end

function get_half_DL_charge(cell::AYALG1iBoltzmann, side="l")
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

    nes = nLAus(S) / nLAu * exp(-data.Ge / kB / data.T) * BoltzmannAu_ne(data, V) # [# electrons/ ISR area]
    
    ISRcontribution = (data.boundary_charge_fac)*e0*(nLAus(S) - nes)
    #@show V, ISRcontribution, -QrAuL[ipsi, bulk]

    return -QrAuL[ipsi, bulk] + ISRcontribution # the sign of the surface charge is opposite, see surface Poisson
end

bulk_charge(cell::AbstractCell) = -1 .*VoronoiFVM.integrate(cell.system, cell.system.physics.reaction, cell.U) # nspec x nregion
boundary_charge(cell::AbstractCell) = -1 .*VoronoiFVM.integrate(cell.system, cell.system.physics.breaction, cell.U, boundary=true)


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