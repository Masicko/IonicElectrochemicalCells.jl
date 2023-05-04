#=
Domain annotation
=#
const Ω_Aul, Ω_YSZ, Ω_Aur = 66, 2, 3
const bulk_domains = (Ω_Aul, Ω_YSZ, Ω_Aur)
const Γ_Aul, Γ_YSZl, Γ_YSZr, Γ_Aur = 1, 3, 4, 2
const interfaces = (Γ_Aul, Γ_YSZl, Γ_YSZr, Γ_Aur)
const ISR = (Γ_YSZl, Γ_YSZr)
const ipsi, iyV = 1, 2
const bulk_species = (ipsi, iyV)
const iyVs = 3
const surface_species = (iyVs)
const species_names = ("ipsi", "iyV")
# extended quantities
const e_Ω_YSZl, e_Ω_YSZr = 2, 4


#=
Grid parameters
=#
const electrolyte_length = 1.0e-8 # YHElectrolyte
const electrolyte_thickness = 1e-3
const electrode_thickness = 1e-5
const dmin = 1.0e-11

#=
YSZ parameters
=#

#
const nu = 0.93 # backward compatibility for the original model

# lattice
const nLYSZ = 7.46268656716418e27 # https://www.msesupplies.com/products/ysz-single-crystal?variant=31250580111418
enLYSZ(S::Float64) = S*nLYSZ

const x_frac = 0.13
const mC = 4.0
const mA = 8.0
# charges and electrostatics
#const epsYSZ = 28.0 * ε0
const epsYSZ = 28.0 * ε0          # https://link.springer.com/article/10.1007/BF01538778
const zV = 2.0
const zi = -zV
""" lattice charge density in YSZ """
zCf(x) = 4.0 * (1.0 - x) / (1.0 + x) + 2.0 * 3.0 * x / (1.0 + x)

const zC = zCf(x_frac)
const zL = mC * zC - 2 * mA

YSZ_charge_density(nV, x) = e0 * nLYSZ * (mC * zCf(x) + mA * zi) + e0 * zV * nV # TODO promote x_frac to YSZ parameters
YSZ_charge_density(nV) = YSZ_charge_density(nV, x_frac) #e0*nLYSZ*(mC*zC + mA*zi + zV*nV/nLYSZ)
electroneutral_nV_YSZ(x) = -(mC * zCf(x) + mA * zi) * nLYSZ / zV
const electroneutral_nV = electroneutral_nV_YSZ(x_frac)
nVmax(alpha::Float64) = nLYSZ*(mA - (zC*mC)/zV*(1-alpha))
nVmax(sys::VoronoiFVM.AbstractSystem) = nVmax(sys.physics.data.alpha)
yV_en(alpha::Float64) = electroneutral_nV / nVmax(alpha)

# extended 
e_YSZ_charge_density(nV, x, S) = e0 * enLYSZ(S) * (mC * zCf(x) + mA * zi) + e0 * zV * nV # TODO promote x_frac to YSZ parameters
e_YSZ_charge_density(nV, S) = e_YSZ_charge_density(nV, x_frac, S)
e_electroneutral_nV_YSZ(x, S) = -(mC * zCf(x) + mA * zi) * enLYSZ(S) / zV
e_electroneutral_nV(S) = e_electroneutral_nV_YSZ(x_frac, S)
e_nVmax(alpha::Float64, S) = enLYSZ(S)*(mA - (zC*mC)/zV*(1-alpha))
e_yV_en(alpha::Float64) = e_electroneutral_nV(S) / e_nVmax(alpha, S)


#=
Au (gold) domain parameters
=#
const epsAu = 6.9 * ε0 # Separation of the contribution of free and bound electrons into real and imaginary parts of the dielectric constant of gold.     Shklyarevskii, I. N.; Pakhmov, P. L. USSR. Optika i Spektroskopiya  (1973),  34(1),  163-6.
const ze = -1.0

# lattice density
const nLAu = 1 / (4 * pi / 3 * (3.01 * a0)^3) # RM Martin, Electronic structure (eq 5.1, )
enLAu(S) = S*nLAu

# charge
Au_charge_density(ne) = e0 * (nLAu * 1.0 + ze * ne)
e_Au_charge_density(ne, S) = e0 * (enLAu(S) * 1.0 + ze * ne)
# electroneutral concentrations
const neutral_yeAU = 1.0
const electroneutral_nAu = neutral_yeAU * nLAu
e_electroneutral_nAu(S) = neutral_yeAU * nLAu(S)

#=
Interface specific reagion (ISR) 
=#


nLYSZs(S::Float64) = (enLYSZ(S))^(2 / 3)  # [# YSZ faces/m^2]
nLAus(S::Float64) = (enLAu(S))^(2 / 3)    # [# Au faces/m^2]

const mCs = 2               # [# cation sites/YSZ face]
const mAs = 4               # [# anion sites/YSZ face]
const zLYSZs = mCs * zC - 2 * mAs # [static elementary charge/YSZ face] <=> fully occupied oxide vacancies

YSZ_surface_charge(yVs, alphas, S::Float64) = e0 * (nVmaxs(alphas, S) * yVs * zV + nLYSZs(S) * zLYSZs)
Au_surface_charge(yes, S::Float64) = e0 * (nLAus(S)*( 1 - yes))    
ISR_staticcharge(S::Float64) = nLYSZs(S) * zLYSZs + 1 * nLAus(S)                          # per ISR area [C/m^2]
ISR_chargedensity(nes, nVs, S::Float64) = e0 * (-nes + nVs * zV + ISR_staticcharge(S)) # per cell cross-section [C/m^2]

# electroneutral with respect to only YSZ surface lattice (electrons are neglected now)
electroneutral_nVs_YSZ(x, S::Float64) = (-mCs * zCf(x)/zV + mAs) * nLYSZs(S)

nVmaxs(a::Float64, S::Float64)::Float64 = nLYSZs(S)*(-zC*mCs / zV * (1 - a) + mAs) # = nVmaxs / S_l nLYSZs(S) = (mAs - (1-a)*mCs*zC/zV)
yVs_en(a::Float64, S::Float64) = electroneutral_nVs_YSZ(x_frac, S)/ nVmaxs(a, S)

ISR_arearatio(bnode,data) = (bnode.region == Γ_YSZl ? data.SL : data.SR)
reduced_voltage(u, bnode, data) = u[ipsi] - (bnode.region == Γ_YSZl ? data.bias : 0.0)
ISR_electrondensity(u, bnode, data) = nLAus(ISR_arearatio(bnode,data)) / nLAu * exp(-data.Ge / kB / data.T) * BoltzmannAu_ne(data, reduced_voltage(u,bnode,data)) # [# electrons/ ISR area]




# model parameters
@kwdef struct SharedParams{A,B}
    T::A = Param(T0 + 800, bounds=(600, 900), name = "temperature")
    bias::B = Param(0.0, bounds=(-3.0, 3.0))
end
@kwdef struct YSZparams{A,B,C}
    AYSZ::A = Param(0.0 * 1.0e-3, bounds=(0.0, 1e2))
    DYSZ::B = Param(1.0e-6, bounds=(1e-5, 1e2))
    alpha::C = Param(0.005, bounds=(0.0, 1.0))
end
@kwdef struct Auparams{A}
    DAu::A = Param(1.0e-8, bounds=(1e-5, 1e2))
end

@kwdef struct ISRparameters{A,B,C,D,E,F,G,H,I,J}
    alphas::A = Param(0.005, bounds=(0.0,1.0)) # [1] ratio of admissible vacancies at ISR
    AYSZs::B = Param(0.0, bounds=(0.0,1.0)) # [eV] interaction energy of vacancies at ISR
    GA::C = Param(0.0 * e0, bounds=(0.0,1.0)) # [eV] Gibbs energy of vacancy adsorption
    Ge::D = Param(0.0 * e0, bounds=(0.0,1.0)) # [eV] Gibbs energy of electron adsorption
    kA::E = Param(1.0e25, bounds=(0.0,1.0)) # [1/m^2/s] rate of vacancy adsorption
    dL::F = Param(1.0 * nm, bounds=(0.0,1.0)) # [nm] thickness of the left ISR
    dR::G = Param(1.0 * nm, bounds=(0.0,1.0)) # [nm] thickness of the right ISR
    SL::H = Param(1.0, bounds=(0.0,1.0)) # [1] (ISR area)/(cross section area) -- left
    SR::I = Param(1.0, bounds=(0.0,1.0)) # [1] (ISR area)/(cross section area) -- right
    boundary_charge_fac::J = Param(1.0, bounds=(-100, 100)) # olny for debugging ... to effectively switch off boundary charge
end