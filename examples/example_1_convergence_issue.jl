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
        :SL => 0.1,
        :SR => 0.7, 
        :T => 800,
        :Ge => 0.0e0,
        :DYSZ => 1.0e-6,
        :kA => 1.0e25
    ]);

    slow_prms = Dict(
        :AYSZ => 1.8279200439453123, 
        :AYSZs => 7.838887329101563, 
        :alpha => 1.913371568830365e-6, 
        :alphas => 6.79157052489046e-7, 
        :GA => 4.635985709513672e-20,
        :SL => 0.49561767578125004, 
        :SR => 0.26644287109375,
        :T => 800.0, 
        :Ge => 0.0, 
        :DYSZ => 1.0e-6, 
        :kA => 1.0e25,  
    );

    slower_prms = Dict(:alpha => 8.228946490285324e-6, :SR => 0.73885498046875, :kA => 1.0e25, :SL => 0.11461181640625001, :AYSZs => 2.694999145507813, :Ge => 0.0, :alphas => 2.7724191290627477e-5, :AYSZ => 4.472120727539063, :bias => 0.0, :T => 800.0, :DYSZ => 1.0e-6, :GA => 6.18183582122461e-20)

    extra_slow_prms =  Dict(:alpha => 1.6175670983146168e-5, :SR => 0.6130615234375, :kA => 1.0e25, :SL => 0.3783935546875, :AYSZs => 5.697553466796875, :Ge => 0.0, :alphas => 6.275261610889185e-6, :AYSZ => 2.037854736328125, :bias => 0.0, :T => 800.0, :DYSZ => 1.0e-6, :GA => 1.5927888802851564e-20)

    #phys_prms = Dict(:alpha => 1.1086636403875841e-5, :SR => 0.2052490234375, :kA => 1.0e25, :SL => 0.36433105468750004, :AYSZs => 3.0728815917968753, :Ge => 0.0, :alphas => 1.2840127491794392e-6, :AYSZ => 1.162964111328125, :bias => 0.0, :T => 800.0, :DYSZ => 1.0e-6, :GA => 5.798502544535157e-20)

    #phys_prms = Dict(:alpha => 1.1086636403875841e-5, :SR => 0.2052490234375, :kA => 1.0e25, :SL => 0.36433105468750004, :AYSZs => 3.0728815917968753, :Ge => 0.0, :alphas => 1.2840127491794392e-6, :AYSZ => 1.162964111328125, :bias => 0.0, :T => 800.0, :DYSZ => 1.0e-6, :GA => 5.798502544535157e-20)
    
    phys_prms = Dict(:alpha => 3.09109465873862e-5, :SR => 0.7818115234375, :kA => 1.0e25, :SL => 0.43464355468750004, :AYSZs => 6.197490966796875, :Ge => 0.0, :alphas => 1.5537524839148867e-5, :AYSZ => 6.537292236328125, :bias => 0.0, :T => 800.0, :DYSZ => 1.0e-6, :GA => 8.802583733285157e-20)

    
    nonconverg_prms = deepcopy(working_prms)
    nonconverg_prms[:alpha] = 7.0e-7

    for prms_sym in [
            #working_prms,
            #nonconverg_prms,
            #slow_prms,
            #extra_slow_prms,
            phys_prms
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