module example_1_convergence_issue

using IonicElectrochemicalCells
using Plots

iec = IonicElectrochemicalCells

function main(;plot_bool = false)
    
    working_prms = Dict([
        :AYSZ => 1.0,
        :AYSZs => 0.1,
        :alpha => 2.0e-5,
        :alphas => 2.0e-5,
        :GA => 0.4*e0, 
        :SR => 0.7, 
        :SL => 0.1,
        :T => 800,
        :Ge => 0.0e0,
        :DYSZ => 1.0e-6,
        :kA => 1.0e25
    ]);

    nonconverg_prms = deepcopy(working_prms)
    nonconverg_prms[:alpha] = 7.0e-7

    for prms_sym in [
            working_prms,
            nonconverg_prms,
            ]
        Testcell = iec.AYA_Sl()
        try 
            stationary_update!(Testcell, prms_sym, tend=1.0e-2)
            testsim_df = iec.capacitance_measurement!(Testcell, bend=0.3, bstep=0.05, tend=1.0e-0)    
            plot_bool && Plots.plot(testsim_df[:, :bias], 1.0 .*testsim_df[:, :capacitance])
            println("ok")
        catch e
            println(e)
        end
    end
end

end # of module