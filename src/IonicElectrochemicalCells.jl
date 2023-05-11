module IonicElectrochemicalCells

using Base: @kwdef
using CompositeStructs
using DataFrames            # TODO bump to examples
using DocStringExtensions   
using ExtendableGrids
using Markdown
using ModelParameters
using VoronoiFVM            # https://github.com/j-fu/VoronoiFVM.jl
using Plots
using Optim

include("auxiliary.jl");            export set_parameters!, ishless, row2dict
include("cell.jl");                 export AbstractCell, CommonCell #, inival, stationary_update!, update_parameters!
include("grids.jl");                export half_electrolyte_1D, full_electrolyte_1D, cell1D, simplefullcell1D
include("finite_volume_fluxes.jl");   export LGS_flux, Fick_flux, PB_flux
include("units.jl");                export kB, Planck_constant, mₑ, e0, ε0, N_Avo, T0, a0, nm

# cells
### Au_YSZ cell
include("./cells/Au_and_YSZ/definitions.jl")
include("./cells/Au_and_YSZ/solver_control.jl")
include("./cells/Au_and_YSZ/AuYSZAu.jl")
include("./cells/Au_and_YSZ/common_functions.jl")
include("./cells/Au_and_YSZ/macro_stuff.jl")
#
include("./cells/Au_and_YSZ/analytical_solution.jl")
#
#include("./cells/Au_and_YSZ/inspecting_tools.jl")

export YHElectrolyte, AYABoltzmann, AYALGBoltzmann, AYALG1iBoltzmann, AYA_Sl
# methods
export inival, biasshow!, biastest!, stationary_update!, update_parameters!, AuL_charge, get_stored_charge

end
