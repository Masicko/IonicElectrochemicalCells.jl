# stable configuration of the solver
const testing = true
const control = VoronoiFVM.NewtonControl()
#control.tol_absolute = 1e-4
control.tol_relative = 1e-8
control.max_iterations = 500
# control.damp_initial = 1e-6
# control.damp_growth = 1.15
# control.verbose=true
# control.force_first_step = true
# control.handle_exceptions= true
control.Δt = 1.0e-4
control.Δt_min = 1.0e-50
control.Δt_max = 1.0
control.Δu_opt = 0.1 # smaller if imprecise
control.verbose = false #testing ? true : false
control.in_memory = false
control.store_all = false

control.tol_round=1.0e-10
control.max_round=4