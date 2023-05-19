module example_1_convergence_issue

using IonicElectrochemicalCells
using Plots

iec = IonicElectrochemicalCells

function main(;plot_bool = false)
    
    middle_peak_prms = Dict([
        :AYSZ => 1.0,
        :AYSZs => 0.1,
        :alpha => 2.0e-5,
        :alphas => 2.0e-5,
        :GA => 0.4*e0, 
        :SL => 0.1,
        :SR => 0.7, 
        :T => 800.0,
        :Ge => 0.0e0,
        :DYSZ => 1.0e-6,
        :kA => 1.0e25
    ]);
    low_error_prms = Dict(
        :AYSZ => 7.8279200439453123, 
        :AYSZs => 7.838887329101563, 
        :alpha => 1.913371568830365e-7, 
        :alphas => 6.79157052489046e-7, 
        :GA => 4.635985709513672e-20,
        :SL => 0.49561767578125004, 
        :SR => 0.26644287109375,
        :T => 800.0, 
        :Ge => 0.0, 
        :DYSZ => 1.0e-6, 
        :kA => 1.0e25,  
    );
    end_peak_prms =  Dict(:alpha => 1.6175670983146168e-5, :SR => 0.6130615234375, :kA => 1.0e25, :SL => 0.3783935546875, :AYSZs => 5.697553466796875, :Ge => 0.0, :alphas => 6.275261610889185e-6, :AYSZ => 2.037854736328125, :bias => 0.0, :T => 800.0, :DYSZ => 1.0e-6, :GA => 1.5927888802851564e-20)
    #nonconv_prms = Dict(:alpha => 5.37309611489754e-5, :SR => 0.164019775390625, :kA => 1.0e25, :SL => 1.454656982421875, :AYSZs => 0.7499298095703124, :Ge => 0.0, :alphas => 1.4098403040439906e-6, :AYSZ => 7.995117797851563, :bias => 0.0, :T => 1073.0, :DYSZ => 1.0e-6, :GA => 4.667278221896485e-20)
    
    
    
    ok_prms = deepcopy(middle_peak_prms)
    ok_prms[:alpha] = 7.0e-7

    for prms_sym in [
            middle_peak_prms,   # works with control.Δt = 1.0e-4 as the only one
            ok_prms,            # works with control.Δt = 1.0e-2
            low_error_prms,     # needs  control.handle_exceptions= true
            end_peak_prms,      # needs  control.handle_exceptions= true
            #nonconv_prms       # nonconvergent at all ... not physically
            ]
        Testcell = iec.AYA_Sl()
        time = @elapsed try
            stationary_update!(Testcell, prms_sym, tend=1.0e-2)
            testsim_df = iec.capacitance_measurement!(Testcell, bend=0.3, bstep=0.05, tend=1.0e-0, testing=true)
            #
            plot_bool && Plots.plot(testsim_df[:, :bias], 1.0 .*testsim_df[:, :capacitance])
            @show testsim_df
            print("ok    ")
        catch e
            print(e, "   ")
        end
        println(" TIME: ",time)
    end
end

end # of module